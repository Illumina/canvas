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
        private double[] gc;
        private double[] sqrtCoverage;

        /*
            // TODO: pass in x and y instead and not pass in the manifest
            IEnumerable<GenomicBin> onTargetBins = manifest == null ? bins : EnrichmentUtilities.GetOnTargetBins(bins, manifest);
            List<double> x = new List<double>();
            List<double> y = new List<double>();
            foreach (var bin in onTargetBins)
            {
                x.Add(bin.GC);
                y.Add(Math.Sqrt(bin.Count)); // Variance stablization
            }
            double medianY = Utilities.Median(y);
            double[] xArr = x.ToArray();
            double[] yArr = y.ToArray();
         */

        /// <summary>
        /// Based on Donavan Cheng's R script (loessnormalize_dev_v2.R)
        /// </summary>
        /// <param name="bandwidth"></param>
        /// <param name="gc">GC content</param>
        /// <param name="coverage">coverage after variance stabilization. Assumed to be normally distributed</param>
        /// <returns></returns>
        private static double objective(double bandwidth, double[] gc, double[] coverage)
        {
            double medianY = Utilities.Median(coverage);

            LoessInterpolator loess = new LoessInterpolator(bandwidth, 0);
            // LOESS
            double[] fittedY = loess.Train(gc, coverage, 1).Fitted;
            double[] normalized = new double[fittedY.Length];
            for (int i = 0; i < normalized.Length; i++)
            {
                normalized[i] = coverage[i] - fittedY[i] + medianY;
            }
            // another LOESS
            return Utilities.StandardDeviation(loess.Train(gc, normalized, 1).Fitted);
        }
    }
}
