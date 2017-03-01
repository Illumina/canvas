using System;
using System.Collections.Generic;
using System.Linq;
using CanvasCommon;
using System.IO;
using Illumina.Common;
using Isas.Manifests.NexteraManifest;
using Isas.SequencingFiles;

namespace CanvasClean
{
    class CanvasClean
    {
        private static readonly int defaultMinNumberOfBinsPerGC = 100;
        private static int minNumberOfBinsPerGCForWeightedMedian = 100;

        /// <summary>
        /// Displays help at the command line.
        /// </summary>
        /// <param name="p">NDesk OptionSet containing command line parameters.</param>
        static void ShowHelp(OptionSet p)
        {
            Console.WriteLine("Usage: CanvasClean.exe [OPTIONS]+");
            Console.WriteLine("Correct bin counts based on genomic parameters");
            Console.WriteLine();
            Console.WriteLine("Options:");
            p.WriteOptionDescriptions(Console.Out);
        }

        /// <summary>
        /// Debugging option - save off a table of counts by GC bin.  
        /// </summary>
        static void DebugPrintCountsByGC(List<SampleGenomicBin> bins, string filePath)
        {
            int[][] HistogramByGC = new int[EnrichmentUtilities.numberOfGCbins][];
            for (int GC = 0; GC < HistogramByGC.Length; GC++) HistogramByGC[GC] = new int[1024];
            foreach (SampleGenomicBin bin in bins)
            {
                if (bin.Count < 0 || bin.Count >= 1024) continue;
                HistogramByGC[bin.GenomicBin.GC][(int)bin.Count]++;
            }
            using (StreamWriter writer = new StreamWriter(filePath))
            {
                // Header line:
                writer.Write("#Bin\\GC");
                for (int GC = 0; GC < HistogramByGC.Length; GC++)
                    writer.Write("GC{0}\t", GC);
                writer.WriteLine();
                // Body:
                for (int count = 0; count < 1024; count++)
                {
                    writer.Write("{0}\t", count);
                    for (int GC = 0; GC < HistogramByGC.Length; GC++)
                        writer.Write("{0}\t", HistogramByGC[GC][count]);
                    writer.WriteLine();
                }
            }
            Console.WriteLine("Wrote counts-by-GC histogram to {0}", filePath);
        }

        /// <summary>
        /// Perform variance stabilization by GC bins.
        /// </summary>
        /// <param name="bins">Bins whose counts are to be normalized.</param>
        static bool NormalizeVarianceByGC(List<SampleGenomicBin> bins, NexteraManifest manifest = null)
        {
            // DebugPrintCountsByGC(bins, "CountsByGCVariance-Before.txt");
            // An array of lists. Each array element (0-100) will hold a list of counts whose bins have the same GC content.
            List<float>[] countsByGC;
            // Will hold all of the autosomal counts present in 'bins'
            List<float> counts;
            EnrichmentUtilities.GetCountsByGC(bins, manifest, out countsByGC, out counts);

            // Estimate quartiles of all bins genomewide
            var globalQuartiles = Utilities.Quartiles(counts);
            // Will hold interquartile range (IQR) separately for each GC bin
            List<float> localIQR = new List<float>(countsByGC.Length);
            // Will hold quartiles separately for each GC bin
            List<Tuple<float, float, float>> localQuartiles = new List<Tuple<float, float, float>>(countsByGC.Length);

            // calculate interquartile range (IQR) for GC bins and populate localQuartiles list
            for (int i = 0; i < countsByGC.Length; i++)
            {
                if (countsByGC[i].Count == 0)
                {
                    localIQR.Add(-1f);
                    localQuartiles.Add(new Tuple<float, float, float>(-1f, -1f, -1f));
                }
                else if (countsByGC[i].Count >= defaultMinNumberOfBinsPerGC)
                {
                    localQuartiles.Add(Utilities.Quartiles(countsByGC[i]));
                    localIQR.Add(localQuartiles[i].Item3 - localQuartiles[i].Item1);
                }
                else
                {
                    List<Tuple<float, float>> weightedCounts = GetWeightedCounts(countsByGC, i);
                    double[] quartiles = Utilities.WeightedQuantiles(weightedCounts, new List<float>() { 0.25f, 0.5f, 0.75f });
                    localQuartiles.Add(new Tuple<float, float, float>((float)quartiles[0], (float)quartiles[1], (float)quartiles[2]));
                    localIQR.Add((float)(quartiles[2] - quartiles[0]));
                }
            }

            // Identify if particular GC bins have IQR twice as large as IQR genomewide 
            float globalIQR = globalQuartiles.Item3 - globalQuartiles.Item1;
            // Holder for GC bins with large IQR (compared to genomewide IQR)
            int significantIQRcounter = 0;
            for (int i = 10; i < 90; i++)
            {
                if (globalIQR < localIQR[i] * 2f)
                    significantIQRcounter++;
            }

            if (significantIQRcounter <= 0)
                return false;

            // Divide each count by the median count of bins with the same GC content
            foreach (SampleGenomicBin bin in bins)
            {
                var scaledLocalIqr = localIQR[bin.GenomicBin.GC] * 0.8f;
                if (globalIQR >= scaledLocalIqr) continue;

                // ratio of GC bins and global IQRs
                float iqrRatio = scaledLocalIqr / globalIQR;
                var medianGCCount = localQuartiles[bin.GenomicBin.GC].Item2;
                bin.Count = medianGCCount + (bin.Count - medianGCCount) / iqrRatio;
            }

            // DebugPrintCountsByGC(bins, "CountsByGCVariance-After.txt");
            return true;
        }

        /// <summary>
        /// In case there are not enough read counts for a GC bin, construct a weighted list of read counts.
        /// Values from target bin i get full weight, the two neighboring bins get half weight, two-away
        /// neighbors get 1/4 weight, etc.
        /// </summary>
        /// <param name="countsByGC"></param>
        /// <param name="gcBin"></param>
        /// <returns></returns>
        private static List<Tuple<float, float>> GetWeightedCounts(List<float>[] countsByGC, int gcBin)
        {
            List<Tuple<float, float>> weightedCounts = new List<Tuple<float, float>>();
            int radius = 0;
            float weight = 1;
            while (weightedCounts.Count < defaultMinNumberOfBinsPerGC)
            {
                int gcWindowEnd = gcBin + radius;
                int gcWindowStart = gcBin - radius;
                if (gcWindowEnd >= countsByGC.Length && gcWindowStart < 0) { break; }

                if (gcWindowEnd < countsByGC.Length)
                {
                    weightedCounts.AddRange(countsByGC[gcWindowEnd].Select(c => Tuple.Create(c, weight)));
                }

                if (gcWindowStart != gcWindowEnd && gcWindowStart >= 0)
                {
                    weightedCounts.AddRange(countsByGC[gcWindowStart].Select(c => Tuple.Create(c, weight)));
                }
                radius++;
                weight /= 2;
            }

            return weightedCounts;
        }

        /// <summary>
        /// Perform GC normalization depending on the mode
        /// </summary>
        /// <param name="bins">Bins whose counts are to be normalized</param>
        /// <param name="manifest"></param>
        /// <param name="mode">GC normalization mode</param>
        static void NormalizeByGC(List<SampleGenomicBin> bins, NexteraManifest manifest, CanvasGCNormalizationMode mode)
        {
            switch (mode)
            {
                case CanvasGCNormalizationMode.MedianByGC:
                    NormalizeByGC(bins, manifest: manifest);
                    break;
                case CanvasGCNormalizationMode.LOESS:
                    var normalizer = new LoessGCNormalizer(bins, manifest, robustnessIter: 0,
                        countTransformer: x => (double)Math.Log(x),
                        invCountTransformer: x => (float)Math.Exp(x));
                    normalizer.Normalize();
                    break;
                default:
                    throw new ApplicationException("Unsupported Canvas GC normalization mode: " + mode.ToString());
            }
        }

        /// <summary>
        /// Perform a simple GC normalization.
        /// </summary>
        /// <param name="bins">Bins whose counts are to be normalized.</param>
        /// <param name="manifest"></param>
        static void NormalizeByGC(List<SampleGenomicBin> bins, NexteraManifest manifest = null)
        {
            // DebugPrintCountsByGC(bins, "CountsByGC-Before.txt");
            // An array of lists. Each array element (0-100) will hold a list of counts whose bins have the same GC content.
            List<float>[] countsByGC;

            // Will hold all of the autosomal counts present in 'bins'
            List<float> counts;
            EnrichmentUtilities.GetCountsByGC(bins, manifest, out countsByGC, out counts);

            double globalMedian = Utilities.Median(counts);
            double?[] medians = new double?[countsByGC.Length];

            // Compute the median count for each GC bin
            for (int gcBinIndex = 0; gcBinIndex < countsByGC.Length; gcBinIndex++)
            {
                if (countsByGC[gcBinIndex].Count >= defaultMinNumberOfBinsPerGC)
                {
                    medians[gcBinIndex] = Utilities.Median(countsByGC[gcBinIndex]);
                }
                else
                {
                    List<Tuple<float, float>> weightedCounts = GetWeightedCounts(countsByGC, gcBinIndex);
                    medians[gcBinIndex] = Utilities.WeightedMedian(weightedCounts);
                }
            }

            // Divide each count by the median count of bins with the same GC content
            for (int gcBinIndex = 0; gcBinIndex < bins.Count; gcBinIndex++)
            {
                double? median = medians[bins[gcBinIndex].GenomicBin.GC];
                if (median != null && median > 0)
                    bins[gcBinIndex].Count = (float)(globalMedian * (double)bins[gcBinIndex].Count / median);
            }
            // DebugPrintCountsByGC(bins, "CountsByGC-After.txt");
        }

        /// <summary>
        /// Remove bins with extreme GC content.
        /// </summary>
        /// <param name="bins">Genomic bins in from which we filter out GC content outliers.</param>
        /// <param name="threshold">Minimum number of bins with the same GC content required to keep a bin.</param>
        /// 
        /// The rationale of this function is that a GC normalization is performed by computing the median count
        /// for each possible GC value. If that count is small, then the corresponding normalization constant
        /// is unstable and we shouldn't use these data.
        static List<SampleGenomicBin> RemoveBinsWithExtremeGC(List<SampleGenomicBin> bins, int threshold, NexteraManifest manifest = null)
        {
            // Will hold outlier-removed bins.
            List<SampleGenomicBin> stripped = new List<SampleGenomicBin>();

            // used to count the number of bins with each possible GC content (0-100)
            int[] counts = new int[EnrichmentUtilities.numberOfGCbins];
            double totalCount = 0;
            foreach (SampleGenomicBin bin in manifest == null ? bins : EnrichmentUtilities.GetOnTargetBins(bins, manifest))
            {

                // We only count autosomal bins because these are the ones we computed normalization factor upon.
                if (!GenomeMetadata.SequenceMetadata.IsAutosome(bin.GenomicBin.Chromosome))
                    continue;

                counts[bin.GenomicBin.GC]++;
                totalCount++;
            }

            int averageCountPerGC = Math.Max(minNumberOfBinsPerGCForWeightedMedian, (int)(totalCount / counts.Length));
            threshold = Math.Min(threshold, averageCountPerGC);
            foreach (SampleGenomicBin bin in bins)
            {
                // Remove outlier (not a lot of bins with the same GC content)
                if (counts[bin.GenomicBin.GC] < threshold)
                    continue;
                stripped.Add(bin);
            }

            return stripped;
        }


        /// <summary>
        /// Calculates Standard Deviation separately for each chromosome and output their average 
        /// </summary>
        public static double GetLocalStandardDeviationAverage(List<double> list, List<string> chromosome)
        {
            List<double> medianAbsoluteDeviations = new List<double>();
            int iStart = 0;
            for (int i = 0; i < list.Count; i++)
            {
                if (chromosome[i] != chromosome[iStart])
                {
                    int iEnd = i; // 0-based, exclusive
                    medianAbsoluteDeviations.Add(Utilities.Mad(list, iStart, iEnd));
                    iStart = i;
                }
            }
            medianAbsoluteDeviations.Add(Utilities.Mad(list, iStart, list.Count));
            return medianAbsoluteDeviations.Average();
        }

        /// <summary>
        /// Estimate local standard deviation (SD).
        /// </summary>
        /// <param name="bins">Genomic bins from which we filter out local SD outliers associated with FFPE biases.</param>
        /// <param name="threshold">Median SD value which is used to determine whereas to run RemoveBinsWithExtremeLocalMad on a sample and which set of bins to remove (set as threshold*5).</param>
        /// The rationale of this function is that standard deviation of difference of consecutive bins values, when taken over a small range of bin (i.e. 20 bins),
        /// has a distinct distribution for FFPE compared to Fresh Frozen (FF) samples. This property is used to flag and remove such bins.

        static double getLocalStandardDeviation(List<SampleGenomicBin> bins)
        {
            // Will hold consecutive bin count difference (approximates Skellam Distribution: mean centred on zero so agnostic to CN changes)
            double[] countsDiffs = new double[bins.Count - 1];

            for (int binIndex = 0; binIndex < bins.Count - 1; binIndex++)
            {
                countsDiffs[binIndex] = Convert.ToDouble(bins[binIndex + 1].Count - bins[binIndex].Count);
            }

            // holder of local SD values (SDs of 20 bins)
            List<double> localSDs = new List<double>();
            List<string> chromosomeBin = new List<string>();

            // calculate local SD metric
            int windowSize = 20;
            for (int windowEnd = windowSize, windowStart = 0; windowEnd < countsDiffs.Length; windowStart += windowSize, windowEnd += windowSize)
            {
                double localSD = Utilities.StandardDeviation(countsDiffs, windowStart, windowEnd);
                localSDs.Add(localSD);
                chromosomeBin.Add(bins[windowStart].GenomicBin.Chromosome);
                for (int binIndex = windowStart; binIndex < windowEnd; binIndex += 1)
                {
                    bins[binIndex].CountDeviation = localSD;
                }
            }

            // average of local SD metric
            double localSDaverage = GetLocalStandardDeviationAverage(localSDs, chromosomeBin);
            return localSDaverage;
        }

        /// <summary>
        /// Remove bin regions with extreme local standard deviation (SD).
        /// Assume that MadOfDiffs has been set in GetLocalStandardDeviationAverage().
        /// </summary>
        /// <param name="bins">Genomic bins from which we filter out local SD outliers associated with FFPE biases.</param>
        /// <param name="threshold">Median SD value which is used to determine whereas to run RemoveBinsWithExtremeLocalMad on a sample and which set of bins to remove (set as threshold*5).</param>
        /// The rationale of this function is that standard deviation of difference of consecutive bins values, when taken over a small range of bin (i.e. 20 bins),
        /// has a distinct distribution for FFPE compared to Fresh Frozen (FF) samples. This property is used to flag and remove such bins.
        static List<SampleGenomicBin> RemoveBinsWithExtremeLocalSD(List<SampleGenomicBin> bins, double localSDaverage, double threshold, string outFile)
        {
            // Will hold FFPE outlier-removed bins 
            List<SampleGenomicBin> strippedBins = new List<SampleGenomicBin>();

            // remove bins with extreme local SD (populating new list is faster than removing from existing one)
            foreach (SampleGenomicBin bin in bins)
            {
                // do not strip bins for samples with local SD metric average less then the threshold
                if (bin.CountDeviation > threshold * 2.0 && localSDaverage > 5.0)
                    continue;
                strippedBins.Add(bin);
            }
            return strippedBins;
        }

        /// <summary>
        /// Removes bins that are genomically large. Typically centromeres and other nasty regions.
        /// </summary>
        /// <param name="bins">Genomic bins.</param>
        static List<SampleGenomicBin> RemoveBigBins(List<SampleGenomicBin> bins)
        {
            List<int> sizes = new List<int>(bins.Count);

            foreach (SampleGenomicBin bin in bins)
                sizes.Add(bin.Size);

            sizes.Sort();

            // Get the 98th percentile of bin sizes
            int index = (int)(0.98 * (double)bins.Count);
            if (index >= sizes.Count)
            {
                Console.Error.WriteLine("Warning in CanvasClean: Too few bins to do outlier removal");
                return bins;
            }
            int thresh = sizes[index];

            List<SampleGenomicBin> stripped = new List<SampleGenomicBin>();

            // Remove bins whose size is greater than the 98th percentile
            foreach (SampleGenomicBin bin in bins)
            {
                if (bin.Size <= thresh)
                    stripped.Add(bin);
            }
            return stripped;
        }

        /// <summary>
        /// Determine if two Poisson counts are unlikely to have come from the same distribution.
        /// </summary>
        /// <param name="a">First count to compare.</param>
        /// <param name="b">Second count to compare.</param>
        /// <returns>True if a and b are unlikely to have arisen from the same Poisson distribution.</returns>
        static bool SignificantlyDifferent(float a, float b)
        {

            double mu = ((double)a + (double)b) / 2;

            if (a + b == 0)
                return false;

            // Calculate Chi-Squared statistic
            double da = (double)a - mu;
            double db = (double)b - mu;
            double chi2 = (da * da + db * db) / mu;

            // Is Chi-Squared greater than the 99th percentile of the Chi-Squared distribution with 1 degree of fredom?
            if (chi2 > 6.635)
                return true;

            return false;
        }

        /// <summary>
        /// Removes point outliers from the dataset.
        /// </summary>
        /// <param name="bins">Genomic bins.</param>
        static List<SampleGenomicBin> RemoveOutliers(List<SampleGenomicBin> bins)
        {
            List<SampleGenomicBin> stripped = new List<SampleGenomicBin>();

            // Check each point to see if it is different than both the point to left and the point to the right
            for (int binIndex = 0; binIndex < bins.Count; binIndex++)
            {
                bool hasPreviousBin = binIndex > 0;
                bool hasNextBin = binIndex < bins.Count - 1;
                string currentBinChromosome = bins[binIndex].GenomicBin.Chromosome;
                string previousBinChromosome = hasPreviousBin ? bins[binIndex - 1].GenomicBin.Chromosome : null;
                string nextBinChromosome = hasNextBin ? bins[binIndex + 1].GenomicBin.Chromosome : null;
                // Different chromosome on both sides
                if ((hasPreviousBin && !currentBinChromosome.Equals(previousBinChromosome))
                    && (hasNextBin && !currentBinChromosome.Equals(nextBinChromosome)))
                    continue;
                // Same chromosome on at least on side or it's the only bin
                if ((hasPreviousBin && bins[binIndex].GenomicBin.Chromosome.Equals(previousBinChromosome) && !SignificantlyDifferent(bins[binIndex].Count, bins[binIndex - 1].Count))
                    || (hasNextBin && bins[binIndex].GenomicBin.Chromosome.Equals(nextBinChromosome) && !SignificantlyDifferent(bins[binIndex].Count, bins[binIndex + 1].Count))
                    || (!hasPreviousBin && !hasNextBin))
                {
                    stripped.Add(bins[binIndex]);
                }
            }

            return stripped;
        }

        static int Main(string[] args)
        {
            Utilities.LogCommandLine(args);
            string inFile = null;
            string outFile = null;
            bool doGCnorm = false;
            bool doSizeFilter = false;
            bool doOutlierRemoval = false;
            string ffpeOutliersFile = null;
            string manifestFile = null;
            CanvasGCNormalizationMode gcNormalizationMode = CanvasGCNormalizationMode.MedianByGC;
            string modeDescription = String.Format("gc normalization mode. Available modes: {0}. Default: {1}",
                String.Join(", ", Enum.GetValues(typeof(CanvasGCNormalizationMode)).Cast<CanvasGCNormalizationMode>()),
                gcNormalizationMode);
            bool needHelp = false;

            OptionSet p = new OptionSet()
            {
                { "i|infile=",        "input file - usually generated by CanvasBin",      v => inFile = v },
                { "o|outfile=",       "text file to output containing cleaned bins",      v => outFile = v },
                { "g|gcnorm",         "perform GC normalization",                         v => doGCnorm = v != null },
                { "s|filtsize",       "filter out genomically large bins",                v => doSizeFilter = v != null },
                { "r|outliers",       "filter outlier points",                            v => doOutlierRemoval = v != null },
                { "f|ffpeoutliers=",   "filter regions of FFPE biases",                   v => ffpeOutliersFile = v },
                { "t|manifest=",      "Nextera manifest file",                            v => manifestFile = v },
                { "w|weightedmedian=", "Minimum number of bins per GC required to calculate weighted median", v => minNumberOfBinsPerGCForWeightedMedian = int.Parse(v) },
                { "m|mode=",          modeDescription,                                    v => gcNormalizationMode = Utilities.ParseCanvasGCNormalizationMode(v) },
                { "h|help",           "show this message and exit",                       v => needHelp = v != null },
            };

            List<string> extraArgs = p.Parse(args);

            if (needHelp)
            {
                ShowHelp(p);
                return 0;
            }

            if (inFile == null || outFile == null)
            {
                ShowHelp(p);
                return 0;
            }

            // Does the input file exist?
            if (!File.Exists(inFile))
            {
                Console.WriteLine("CanvasClean.exe: File {0} does not exist! Exiting.", inFile);
                return 1;
            }

            List<SampleGenomicBin> bins = CanvasIO.ReadFromTextFile(inFile);

            if (doOutlierRemoval)
                bins = RemoveOutliers(bins);

            if (doSizeFilter)
                bins = RemoveBigBins(bins);

            // do not run FFPE outlier removal on targeted/low coverage data
            if (ffpeOutliersFile != null && bins.Count < 50000)
            {
                ffpeOutliersFile = null;
            }

            // estimate localSD metric to use in doFFPEOutlierRemoval later and write to a text file 
            double LocalSD = -1.0;
            if (ffpeOutliersFile != null)
            {
                LocalSD = getLocalStandardDeviation(bins);
                CanvasIO.WriteLocalSDToTextFile(ffpeOutliersFile, LocalSD);
            }

            if (doGCnorm)
            {
                NexteraManifest manifest = manifestFile == null ? null : new NexteraManifest(manifestFile, null, Console.WriteLine);
                List<SampleGenomicBin> strippedBins = gcNormalizationMode == CanvasGCNormalizationMode.MedianByGC
                    ? RemoveBinsWithExtremeGC(bins, defaultMinNumberOfBinsPerGC, manifest: manifest)
                    : bins;
                if (strippedBins.Count == 0)
                {
                    Console.Error.WriteLine("Warning in CanvasClean: Coverage too low to perform GC correction; proceeding without GC correction");
                }
                else
                {
                    bins = strippedBins;
                    NormalizeByGC(bins, manifest, gcNormalizationMode);
                    // Use variance normalization only on large exome panels and whole genome sequencing
                    // The treshold is set to 10% of an average number of bins on CanvasClean data
                    if (ffpeOutliersFile != null && bins.Count > 500000)
                    {
                        bool isNormalizeVarianceByGC = NormalizeVarianceByGC(bins, manifest: manifest);
                        // If normalization by variance was run (isNormalizeVarianceByGC), perform mean centering by using NormalizeByGC 
                        if (isNormalizeVarianceByGC)
                            NormalizeByGC(bins, manifest, gcNormalizationMode);
                    }

                }
            }

            if (ffpeOutliersFile != null)
            {
                // threshold 20 is derived to separate FF and noisy FFPE samples (derived from a training set of approx. 40 samples)
                List<SampleGenomicBin> LocalMadstrippedBins = RemoveBinsWithExtremeLocalSD(bins, LocalSD, 20, outFile);
                bins = LocalMadstrippedBins;
            }

            CanvasIO.WriteToTextFile(outFile, bins);
            return 0;

        }

    }

}
