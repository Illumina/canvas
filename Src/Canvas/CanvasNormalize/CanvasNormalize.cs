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
            
            GetWeightedAverageBinCount(parameters.normalBedPaths, parameters.weightedAverageNormalBedPath, manifest: manifest);
            
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
                    List<double> binCounts;
                    // If a manifest is available, get the median of bins overlapping the targeted regions only.
                    // For small panels, there could be a lot of bins with zero count and the median would be 0 if taken over all the bins, resulting in division by zero.
                    double median = manifest == null ? GetMedianBinCount(binnedPath, out binCounts) :
                        GetMedianBinCount(binnedPath, out binCounts, manifest);
                    weights[sampleIndex] = median > 0 ? 1.0 / median : 0;
                    binCountsBySample[sampleIndex] = binCounts;
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

        private static double GetMedianBinCount(string binnedPath, out List<double> binCounts)
        {
            binCounts = new List<double>();
            using (GzipReader reader = new GzipReader(binnedPath))
            {
                string line;
                string[] toks;
                while ((line = reader.ReadLine()) != null)
                {
                    toks = line.Split('\t');
                    binCounts.Add(double.Parse(toks[3]));
                }
            }
            return  Utilities.Median(binCounts); // median bin count of a sample across bins
        }

        private static double GetMedianBinCount(string binnedPath, out List<double> binCounts, NexteraManifest manifest)
        {
            List<double> tmpBinCounts = new List<double>();
            List<int> onTargetIndices = new List<int>();

            var regionsByChrom = manifest.GetManifestRegionsByChromosome();
            string currChrom = null;
            List<NexteraManifest.ManifestRegion> regions = null; // 1-based regions
            int regionIndex = -1;
            bool onTarget = false;
            using (GzipReader reader = new GzipReader(binnedPath))
            {
                string line;
                string[] toks;
                int binIdx = 0;
                while ((line = reader.ReadLine()) != null)
                {
                    toks = line.Split('\t');
                    string chrom = toks[0];
                    int start = int.Parse(toks[1]); // 0-based, inclusive
                    int stop = int.Parse(toks[2]); // 0-based, exclusive
                    if (currChrom != chrom)
                    {
                        currChrom = chrom;
                        onTarget = false;
                        if (!regionsByChrom.ContainsKey(currChrom))
                        {
                            regions = null;
                        }
                        else
                        {
                            regions = regionsByChrom[currChrom];
                            regionIndex = 0;
                        }
                    }
                    while (regions != null && regionIndex < regions.Count && regions[regionIndex].End < start + 1)
                    {
                        regionIndex++;
                    }
                    if (regions != null && regionIndex < regions.Count && regions[regionIndex].Start <= stop) // overlap
                    {
                        onTarget = true;
                    }
                    else
                    {
                        onTarget = false;
                    }

                    if (onTarget) { onTargetIndices.Add(binIdx); }

                    tmpBinCounts.Add(double.Parse(toks[3]));
                    binIdx++;
                }
            }
            binCounts = tmpBinCounts;
            // median bin count of a sample across on-target bins
            return Utilities.Median(onTargetIndices.Select(idx => tmpBinCounts[idx]).ToList());
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
