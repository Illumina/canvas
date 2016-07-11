using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Isas.Shared;
using SequencingFiles;

namespace CanvasNormalize
{
    public class BestLR2ReferenceGenerator : IReferenceGenerator
    {
        private readonly IFileLocation _sampleBinnedFile;
        private readonly IEnumerable<IFileLocation> _controlBinnedFiles;
        private readonly NexteraManifest _manifest;

        public BestLR2ReferenceGenerator(IFileLocation sampleBinnedFile, IEnumerable<IFileLocation> controlBinnedFiles,
            NexteraManifest manifest)
        {
            if (!sampleBinnedFile.Exists)
                throw new FileNotFoundException(sampleBinnedFile.FullName + " does not exist.");
            foreach (var binnedFile in controlBinnedFiles)
                if (!binnedFile.Exists)
                    throw new FileNotFoundException(binnedFile.FullName + " does not exist.");

            _sampleBinnedFile = sampleBinnedFile;
            _controlBinnedFiles = controlBinnedFiles;
            _manifest = manifest;
        }

        public void Run(IFileLocation outputFile)
        {
            int bestNormalSampleIndex = 0;
            int normalSampleCount = _controlBinnedFiles.Count();
            if (normalSampleCount > 1) // find the best normal
            {
                List<double[]> binCountsByNormalSample = new List<double[]>();
                for (int normalSampleIndex = 0; normalSampleIndex < normalSampleCount; normalSampleIndex++)
                {
                    var controlBinnedFile = _controlBinnedFiles.ElementAt(normalSampleIndex);
                    var binCounts = new BinCounts(controlBinnedFile.FullName, manifest: _manifest);
                    List<double> counts = binCounts.OnTargetCounts;
                    double median = binCounts.OnTargetMedianBinCount;
                    // If a manifest is available, get the median of bins overlapping the targeted regions only.
                    // For small panels, there could be a lot of bins with zero count and the median would be 0 if taken over all the bins, resulting in division by zero.
                    double weight = median > 0 ? 1.0 / median : 0;
                    binCountsByNormalSample.Add(counts.Select(cnt => cnt * weight).ToArray());
                }
                double[] tumorBinCounts;
                {
                    var binCounts = new BinCounts(_sampleBinnedFile.FullName, manifest: _manifest);
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
            var srcBinnedFile = _controlBinnedFiles.ElementAt(bestNormalSampleIndex);
            if (outputFile.Exists) { outputFile.Delete(); }
            srcBinnedFile.CopyTo(outputFile);
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
    }
}
