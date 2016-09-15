using System;
using System.Collections.Generic;
using System.Linq;
using CanvasCommon;
using Isas.SequencingFiles;

namespace CanvasClean
{
    public class LoessGCNormalizer
    {
        private double[] gcs;
        private double[] counts;
        private List<int> withoutChrY;
        private IEnumerable<GenomicBin> bins;
        private NexteraManifest manifest;
        private Func<float, double> countTransformer = x => (double)x;
        private Func<double, float> invCountTransformer = x => (float)x;
        private int robustnessIter = 0;


        public LoessGCNormalizer(IEnumerable<GenomicBin> bins, NexteraManifest manifest, int robustnessIter = 2,
            Func<float, double> countTransformer = null, Func<double, float> invCountTransformer = null)
        {
            this.bins = bins;
            this.manifest = manifest;
            if (robustnessIter >= 0) { this.robustnessIter = robustnessIter; }
            if (countTransformer != null && invCountTransformer != null)
            {
                this.countTransformer = countTransformer;
                this.invCountTransformer = invCountTransformer;
            }
            initialize();
        }

        private void initialize()
        {
            IEnumerable<GenomicBin> onTargetBins = manifest == null ? bins : EnrichmentUtilities.GetOnTargetBins(bins, manifest);

            List<double> x = new List<double>();
            List<double> y = new List<double>();
            withoutChrY = new List<int>();
            int i = 0; // index into x and y
            foreach (var bin in onTargetBins)
            {
                double count = countTransformer(bin.CountBin.Count); // Variance stablization
                if (!double.IsInfinity(count))
                {
                    x.Add(bin.GC);
                    y.Add(count);
                    string chrom = bin.Chromosome.ToLower();
                    bool isChrY = chrom == "chry" || chrom == "y";
                    if (!isChrY) { withoutChrY.Add(i); }
                    i++;
                }
            }

            gcs = x.ToArray();
            counts = y.ToArray();
        }

        public void Normalize()
        {
            // Find the best bandwidth without chrY
            double[] gcsNoChrY = withoutChrY.Select(i => gcs[i]).ToArray();
            double[] countsNoChrY = withoutChrY.Select(i => counts[i]).ToArray();
            double bestBandwidth = findBestBandwith(0.3, 0.75, gcsNoChrY, countsNoChrY);

            // Fit LOESS
            double medianY = Utilities.Median(counts);
            int minGC = (int)gcs.Min();
            int maxGC = (int)gcs.Max();
            LoessInterpolator loess = new LoessInterpolator(bestBandwidth, 0);
            var model = loess.Train(gcs, counts, 1, computeFitted: false);
            double[] fittedByGC = model.Predict(Enumerable.Range(minGC, maxGC).Select(i => (double)i));
            // Smooth
            foreach (GenomicBin bin in bins)
            {
                int i = Math.Min(fittedByGC.Length - 1, Math.Max(0, bin.GC - minGC));
                double smoothed = countTransformer(bin.CountBin.Count) - fittedByGC[i] + medianY;
                bin.CountBin.Count = invCountTransformer(smoothed);
            }
        }

        private static double findBestBandwith(double minBandwidth, double maxBandwidth, double[] gcs, double[] counts)
        {
            minBandwidth = Math.Max(2.0 / gcs.Length, minBandwidth);
            maxBandwidth = Math.Min(1.0, maxBandwidth);
            if (maxBandwidth < minBandwidth) { maxBandwidth = minBandwidth; }

            return Utilities.GoldenSectionSearch(b => objective(b, gcs, counts), minBandwidth, maxBandwidth);
        }

        /// <summary>
        /// Based on Donavan Cheng's R script (loessnormalize_dev_v2.R)
        /// </summary>
        /// <param name="bandwidth"></param>
        /// <param name="gc">GC content</param>
        /// <param name="coverage">coverage after variance stabilization. Assumed to be normally distributed</param>
        /// <returns></returns>
        private static double objective(double bandwidth, double[] gcs, double[] counts)
        {
            double medianY = Utilities.Median(counts);
            int minGC = (int)gcs.Min();
            int maxGC = (int)gcs.Max();

            LoessInterpolator loess = new LoessInterpolator(bandwidth, 0);
            // LOESS
            double[] normalized = new double[counts.Length];
            {
                var model = loess.Train(gcs, counts, 1, computeFitted: false);
                double[] fittedByGC = model.Predict(Enumerable.Range(minGC, maxGC).Select(i => (double)i));
                for (int i = 0; i < normalized.Length; i++)
                {
                    int gc = (int)gcs[i];
                    normalized[i] = counts[i] - fittedByGC[gc - minGC] + medianY;
                }
            }
            // another LOESS
            double[] fitted = new double[counts.Length];
            {
                var model = loess.Train(gcs, normalized, 1, computeFitted: false);
                double[] fittedByGC = model.Predict(Enumerable.Range(minGC, maxGC).Select(i => (double)i));
                for (int i = 0; i < fitted.Length; i++)
                {
                    int gc = (int)gcs[i];
                    fitted[i] = fittedByGC[gc - minGC];
                }
            }

            return Utilities.StandardDeviation(fitted);
        }
    }
}
