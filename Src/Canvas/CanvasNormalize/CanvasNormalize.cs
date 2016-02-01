using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using CanvasCommon;
using SequencingFiles;

namespace CanvasNormalize
{
    class CanvasNormalize
    {
        private static readonly double CanvasDiploidBinRatioFactor = 40;

        public static int Run(CanvasNormalizeParameters parameters) 
        {
            NexteraManifest manifest = string.IsNullOrEmpty(parameters.manifestPath) ? null : new NexteraManifest(parameters.manifestPath, null, Console.WriteLine);

            switch (parameters.normalizationMode)
            {
                case CanvasNormalizeMode.BestLR2:
                    GetBestLR2BinCount(parameters.tumorBedPath, parameters.normalBedPaths, parameters.weightedAverageNormalBedPath,
                        manifest: manifest);
                    break;
                case CanvasNormalizeMode.WeightedAverage:
                    GetWeightedAverageBinCount(parameters.normalBedPaths, parameters.weightedAverageNormalBedPath, manifest: manifest);
                    break;
                default:
                    throw new Exception(string.Format("Invalid CanvasNormalize mode '{0}'", parameters.normalizationMode));
            }
            
            GetBinRatio(parameters.tumorBedPath, parameters.weightedAverageNormalBedPath, parameters.outBedPath, parameters.ploidyBedPath);

            return 0;
        }

        private static void GetWeightedAverageBinCount(IEnumerable<string> binnedPaths, string mergedBinnedPath,
            NexteraManifest manifest = null)
        {
            int sampleCount = binnedPaths.Count();
            if (sampleCount == 1) // copy file
            {
                if (File.Exists(binnedPaths.First()))
                {
                    if (File.Exists(mergedBinnedPath)) { File.Delete(mergedBinnedPath); }
                    File.Copy(binnedPaths.First(), mergedBinnedPath);
                }
            }
            else // merge normal samples
            {
                double[] weights = new double[sampleCount];
                List<double>[] binCountsBySample = new List<double>[sampleCount];
                for (int sampleIndex = 0; sampleIndex < sampleCount; sampleIndex++)
                {
                    string binnedPath = binnedPaths.ElementAt(sampleIndex);
                    var binCounts = new BinCounts(binnedPath, manifest: manifest);
                    List<double> counts = binCounts.AllCounts;
                    // If a manifest is available, get the median of bins overlapping the targeted regions only.
                    // For small panels, there could be a lot of bins with zero count and the median would be 0 if taken over all the bins, resulting in division by zero.
                    double median = binCounts.OnTargetMedianBinCount;
                    weights[sampleIndex] = median > 0 ? 1.0 / median : 0;
                    binCountsBySample[sampleIndex] = counts;
                }
                double weightSum = weights.Sum();
                for (int i = 0; i < sampleCount; i++) { weights[i] /= weightSum; } // so weights sum to 1

                // Computed weighted average of bin counts across samples
                using (GzipReader reader = new GzipReader(binnedPaths.First()))
                using (GzipWriter writer = new GzipWriter(mergedBinnedPath))
                {
                    string line;
                    string[] toks;
                    int lineIdx = 0;
                    while ((line = reader.ReadLine()) != null)
                    {
                        toks = line.Split('\t');
                        double weightedBinCount = 0;
                        for (int i = 0; i < sampleCount; i++) { weightedBinCount += weights[i] * binCountsBySample[i][lineIdx]; }
                        toks[3] = String.Format("{0}", weightedBinCount);
                        writer.WriteLine(String.Join("\t", toks));
                        lineIdx++;
                    }
                }
            }
        }

        /// <summary>
        /// Pick the best normal control that has the smallest mean squared log-ratios (LR2s).
        /// </summary>
        /// <param name="tumorBinnedPath"></param>
        /// <param name="normalBinnedPaths"></param>
        /// <param name="bestBinnedPath"></param>
        /// <param name="manifest"></param>
        private static void GetBestLR2BinCount(string tumorBinnedPath, IEnumerable<string> normalBinnedPaths, string bestBinnedPath,
            NexteraManifest manifest = null)
        {
            int bestNormalSampleIndex = 0;
            int normalSampleCount = normalBinnedPaths.Count();
            if (normalSampleCount > 1) // find the best normal
            {
                List<double[]> binCountsByNormalSample = new List<double[]>();
                for (int normalSampleIndex = 0; normalSampleIndex < normalSampleCount; normalSampleIndex++)
                {
                    string normalBinnedPath = normalBinnedPaths.ElementAt(normalSampleIndex);
                    var binCounts = new BinCounts(normalBinnedPath, manifest: manifest);
                    List<double> counts = binCounts.OnTargetCounts;
                    double median = binCounts.OnTargetMedianBinCount;
                    // If a manifest is available, get the median of bins overlapping the targeted regions only.
                    // For small panels, there could be a lot of bins with zero count and the median would be 0 if taken over all the bins, resulting in division by zero.
                    double weight = median > 0 ? 1.0 / median : 0;
                    binCountsByNormalSample.Add(counts.Select(cnt => cnt * weight).ToArray());
                }
                double[] tumorBinCounts;
                {
                    var binCounts = new BinCounts(tumorBinnedPath, manifest: manifest);
                    List<double> counts = binCounts.OnTargetCounts;
                    double tumorMedian = binCounts.OnTargetMedianBinCount;
                    double tumorWeight = tumorMedian > 0 ? 1.0 / tumorMedian : 0;
                    tumorBinCounts = counts.Select(cnt => cnt * tumorWeight).ToArray();
                }

                // Find the best normal sample
                bestNormalSampleIndex = -1;
                double minMeanSquaredLogRatios = double.PositiveInfinity;
                for (int normalSampleIndex = 0; normalSampleIndex < normalSampleCount; normalSampleIndex++)
                {
                    // Get the sum of squared log ratios
                    var result = GetMeanSquaredLogRatios(tumorBinCounts, binCountsByNormalSample[normalSampleIndex]);
                    double meanSquaredLogRatios = result.Item1;
                    int ignoredBinCount = result.Item2;
                    // TODO: Skip a (bad) normal sample if too many bins were ignored.
                    //       Donavan's script skips a normal sample if more than 100 log ratios is NA.
                    //       The cut-off is likely panel-dependent.
                    if (meanSquaredLogRatios < minMeanSquaredLogRatios)
                    {
                        minMeanSquaredLogRatios = meanSquaredLogRatios;
                        bestNormalSampleIndex = normalSampleIndex;
                    }
                }
            }

            // copy file
            string srcBinnedPath = normalBinnedPaths.ElementAt(bestNormalSampleIndex);
            if (File.Exists(srcBinnedPath))
            {
                if (File.Exists(bestBinnedPath)) { File.Delete(bestBinnedPath); }
                File.Copy(srcBinnedPath, bestBinnedPath);
            }
        }

        private static Tuple<double, int> GetMeanSquaredLogRatios(double[] tumorBinCounts, double[] normalBinCounts)
        {
            double sumOfSquaredLogRatios = 0;
            int ignoredBinCount = 0;
            int nBins = 0;
            for (int binIndex = 0; binIndex < tumorBinCounts.Length; binIndex++)
            {
                double tumorBinCount = tumorBinCounts[binIndex];
                double normalBinCount = normalBinCounts[binIndex];

                if (normalBinCount <= 0)
                {
                    ignoredBinCount++;
                    continue;
                }

                double squaredLogRatio;
                try
                {
                    double logRatio = Math.Log(tumorBinCount / normalBinCount);
                    squaredLogRatio = logRatio * logRatio;
                }
                catch (Exception e)
                {
                    Console.WriteLine("Error calculating squared log ratio: {0}. Ignoring bin {1}.", e.Message, binIndex);
                    ignoredBinCount++;
                    continue;
                }

                if (double.IsInfinity(squaredLogRatio) || double.IsNaN(squaredLogRatio))
                {
                    ignoredBinCount++;
                    continue;
                }

                sumOfSquaredLogRatios += squaredLogRatio;
                nBins++;
            }

            double meanSquaredLogRatios = nBins > 0 ? sumOfSquaredLogRatios / nBins : sumOfSquaredLogRatios;
            return Tuple.Create(meanSquaredLogRatios, ignoredBinCount);
        }

        private static void GetBinRatio(string tumorBinnedPath, string normalBinnedPath, string ratioBinnedPath, string ploidyBedPath)
        {
            PloidyInfo referencePloidy = String.IsNullOrEmpty(ploidyBedPath) ? null : PloidyInfo.LoadPloidyFromBedFile(ploidyBedPath);

            using (GzipReader tumorReader = new GzipReader(tumorBinnedPath))
            using (GzipReader normalReader = new GzipReader(normalBinnedPath))
            using (GzipWriter writer = new GzipWriter(ratioBinnedPath))
            {
                string normalLine;
                string tumorLine;
                string[] normalToks;
                string[] tumorToks;
                double normalCount;
                double tumorCount;
                double ratio;
                while ((normalLine = normalReader.ReadLine()) != null)
                {
                    tumorLine = tumorReader.ReadLine();
                    normalToks = normalLine.Split('\t');
                    tumorToks = tumorLine.Split('\t');
                    normalCount = double.Parse(normalToks[3]);
                    tumorCount = double.Parse(tumorToks[3]);
                    // The weighted average count of a bin could be less than 1.
                    // Using these small counts for coverage normalization creates large ratios.
                    // It would be better to just drop these bins so we don't introduce too much noise into segmentation and CNV calling.
                    if (normalCount < 1) { continue; } // skip the bin
                    string chrom = normalToks[0];
                    int start = int.Parse(normalToks[1]);
                    int end = int.Parse(normalToks[2]);
                    // get the normal ploidy from intervalsWithPloidyByChrom
                    double factor = CanvasDiploidBinRatioFactor * GetPloidy(referencePloidy, chrom, start, end) / 2.0;
                    ratio = tumorCount / normalCount * factor;
                    normalToks[3] = String.Format("{0}", ratio);
                    writer.WriteLine(String.Join("\t", normalToks));
                }
            }
        }

        private static int GetPloidy(PloidyInfo referencePloidy, string chrom, int start, int end, int defaultPloidy = 2)
        {
            if (referencePloidy == null) { return defaultPloidy; }

            CanvasSegment segment = new CanvasSegment(chrom, start, end, new List<float>());

            return referencePloidy.GetReferenceCopyNumber(segment);
        }
    }
}
