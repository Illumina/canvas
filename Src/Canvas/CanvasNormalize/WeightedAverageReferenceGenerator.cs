using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Isas.SequencingFiles;
using Isas.Shared.Utilities.FileSystem;

namespace CanvasNormalize
{
    public class WeightedAverageReferenceGenerator : IReferenceGenerator
    {
        private IEnumerable<IFileLocation> _controlBinnedFiles;
        private NexteraManifest _manifest;

        public WeightedAverageReferenceGenerator(IEnumerable<IFileLocation> controlBinnedFiles, NexteraManifest manifest)
        {
            foreach (var binnedFile in controlBinnedFiles)
                if (!binnedFile.Exists)
                    throw new FileNotFoundException(binnedFile.FullName + " does not exist.");

            _controlBinnedFiles = controlBinnedFiles;
            _manifest = manifest;
        }

        public void Run(IFileLocation outputFile)
        {
            int sampleCount = _controlBinnedFiles.Count();
            if (sampleCount == 1) // copy file
            {
                if (outputFile.Exists) { outputFile.Delete(); }
                _controlBinnedFiles.First().CopyTo(outputFile);
            }
            else // merge normal samples
            {
                double[] weights = new double[sampleCount];
                List<double>[] binCountsBySample = new List<double>[sampleCount];
                for (int sampleIndex = 0; sampleIndex < sampleCount; sampleIndex++)
                {
                    var binnedFile = _controlBinnedFiles.ElementAt(sampleIndex);
                    var binCounts = new BinCounts(binnedFile.FullName, manifest: _manifest);
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
                using (GzipReader reader = new GzipReader(_controlBinnedFiles.First().FullName))
                using (GzipWriter writer = new GzipWriter(outputFile.FullName))
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
    }
}
