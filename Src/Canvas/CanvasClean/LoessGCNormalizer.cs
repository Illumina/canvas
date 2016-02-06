using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CanvasCommon;
using SequencingFiles;

namespace CanvasClean
{
    public class LoessGCNormalizer
    {
        private double[] gcs;
        private double[] sqrtCounts;
        private List<int> withoutChrY;
        private IEnumerable<GenomicBin> bins;
        private NexteraManifest manifest;

        public LoessGCNormalizer(IEnumerable<GenomicBin> bins, NexteraManifest manifest)
        {
            this.bins = bins;
            this.manifest = manifest;
            initialize();
        }

        private void initialize()
        {
            IEnumerable<GenomicBin> onTargetBins = manifest == null ? bins : EnrichmentUtilities.GetOnTargetBins(bins, manifest);

            List<double> x = new List<double>();
            List<double> y = new List<double>();
            withoutChrY = new List<int>();
            int i = 0;
            foreach (var bin in onTargetBins)
            {
                x.Add(bin.GC);
                y.Add(Math.Sqrt(bin.Count)); // Variance stablization
                string chrom = bin.Chromosome.ToLower();
                bool isChrY = chrom == "chry" || chrom == "y";
                if (!isChrY) { withoutChrY.Add(i); }
                i++;
            }

            gcs = x.ToArray();
            sqrtCounts = y.ToArray();
        }

        public void Normalize()
        {
            // Find the best bandwidth without chrY
            double[] gcsNoChrY = withoutChrY.Select(i => gcs[i]).ToArray();
            double[] sqrtCountsNoChrY = withoutChrY.Select(i => sqrtCounts[i]).ToArray();
            double bestBandwidth = findBestBandwith(0.3, 0.75, gcsNoChrY, sqrtCountsNoChrY);

            // Fit LOESS
            double medianY = Utilities.Median(sqrtCounts);
            int minGC = (int)gcs.Min();
            int maxGC = (int)gcs.Max();
            LoessInterpolator loess = new LoessInterpolator(bestBandwidth, 0);
            var model = loess.Train(gcs, sqrtCounts, 1, computeFitted: false);
            double[] fittedByGC = model.Predict(Enumerable.Range(minGC, maxGC).Select(i => (double)i));
            // Smooth
            foreach (GenomicBin bin in bins)
            {
                int i = Math.Min(fittedByGC.Length - 1, Math.Max(0, bin.GC - minGC));
                bin.Count = (float)Math.Pow(Math.Sqrt(bin.Count) - fittedByGC[i] + medianY, 2);
            }
        }

        private static double findBestBandwith(double minBandwidth, double maxBandwidth, double[] gcs, double[] sqrtCounts)
        {
            minBandwidth = Math.Max(2.0 / gcs.Length, minBandwidth);
            maxBandwidth = Math.Min(1.0, maxBandwidth);
            if (maxBandwidth < minBandwidth) { maxBandwidth = minBandwidth; }

            return Utilities.GoldenSectionSearch(b => objective(b, gcs, sqrtCounts), minBandwidth, maxBandwidth);
        }

        /// <summary>
        /// Based on Donavan Cheng's R script (loessnormalize_dev_v2.R)
        /// </summary>
        /// <param name="bandwidth"></param>
        /// <param name="gc">GC content</param>
        /// <param name="coverage">coverage after variance stabilization. Assumed to be normally distributed</param>
        /// <returns></returns>
        private static double objective(double bandwidth, double[] gcs, double[] sqrtCounts)
        {
            double medianY = Utilities.Median(sqrtCounts);
            int minGC = (int)gcs.Min();
            int maxGC = (int)gcs.Max();

            LoessInterpolator loess = new LoessInterpolator(bandwidth, 0);
            // LOESS
            double[] normalized = new double[sqrtCounts.Length];
            {
                var model = loess.Train(gcs, sqrtCounts, 1, computeFitted: false);
                double[] fittedByGC = model.Predict(Enumerable.Range(minGC, maxGC).Select(i => (double)i));
                for (int i = 0; i < normalized.Length; i++)
                {
                    int gc = (int)gcs[i];
                    normalized[i] = sqrtCounts[i] - fittedByGC[gc - minGC] + medianY;
                }
            }
            // another LOESS
            double[] fitted = new double[sqrtCounts.Length];
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
