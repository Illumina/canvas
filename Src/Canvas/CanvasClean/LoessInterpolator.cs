using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CanvasCommon;

namespace CanvasClean
{
    public class LoessInterpolator
    {
        public static double DEFAULT_BANDWIDTH = 0.3;
        public static int DEFAULT_ROBUSTNESS_ITERS = 2;

        private readonly int polynomialDegree = 1; // Supports only degree 1 polynomial for now
        /**
         * The bandwidth parameter: when computing the loess fit at
         * a particular point, this fraction of source points closest
         * to the current point is taken into account for computing
         * a least-squares regression.
         * 
         * A sensible value is usually 0.25 to 0.5.
         */
        private double bandwidth;

        /**
         * The number of robustness iterations parameter: this many
         * robustness iterations are done.
         * 
         * A sensible value is usually 0 (just the initial fit without any
         * robustness iterations) to 4.
         */
        private int robustnessIters;

        public LoessInterpolator()
        {
            this.bandwidth = DEFAULT_BANDWIDTH;
            this.robustnessIters = DEFAULT_ROBUSTNESS_ITERS;
        }

        public LoessInterpolator(double bandwidth, int robustnessIters)
        {
            if (bandwidth <= 0 || bandwidth > 1)
            {
                throw new ApplicationException(string.Format("bandwidth must be in the interval (0, 1], but got {0}", bandwidth));
            }
            this.bandwidth = bandwidth;
            if (robustnessIters < 0)
            {
                throw new ApplicationException(string.Format("the number of robustness iterations must be non-negative, but got {0}", robustnessIters));
            }
            this.robustnessIters = robustnessIters;
        }

        /// <summary>
        /// Train a LOESS model
        /// </summary>
        /// <param name="xvals">x values</param>
        /// <param name="yvals">y values</param>
        /// <param name="xStep">x-value step size</param>
        /// <param name="computeFitted">to compute the fitted y values? Always true when robustnessIters > 0</param>
        /// <returns>LOESS model</returns>
        public LoessModel Train(double[] xvals, double[] yvals, double xStep, bool computeFitted = true)
        {
            int[] ascendingOrder = Enumerable.Range(0, xvals.Length).OrderBy(i => xvals[i]).ToArray(); // argsort xvals

            double[] sortedXval = new double[xvals.Length];
            double[] sortedYval = new double[yvals.Length];
            for (int i = 0; i < ascendingOrder.Length; i++)
            {
                sortedXval[i] = xvals[ascendingOrder[i]];
                sortedYval[i] = yvals[ascendingOrder[i]];
            }

            var model = train(sortedXval, sortedYval, xStep, ascendingOrder, computeFitted: computeFitted);

            return model;
        }

        /// <summary>
        /// Train a LOESS model
        /// </summary>
        /// <param name="xval">x values. Assumed sorted in ascending order</param>
        /// <param name="yvals">y values</param>
        /// <param name="xStep">x-value step size</param>
        /// <returns>LOESS model</returns>
        private LoessModel train(double[] xvals, double[] yvals, double xStep, int[] ascendingOrder,
            bool computeFitted = true)
        {
            if (xvals.Length != yvals.Length)
            {
                throw new ApplicationException(string.Format("Loess expects the abscissa and ordinate arrays to be of the same Size, but got {0} abscisssae and {1} ordinatae", xvals.Length, yvals.Length));
            }

            int n = xvals.Length;
            if (n <= 1)
            {
                return new LoessModel(xvals, yvals, new double[] { yvals[0] }, null, ascendingOrder, null);
            }

            checkAllFiniteReal(xvals, true);
            checkAllFiniteReal(yvals, false);
            checkIncreasing(xvals);

            int bandwidthInPoints = (int)Math.Ceiling(bandwidth * n);

            int leastBandwithInPoints = 1 + polynomialDegree;
            if (bandwidthInPoints < leastBandwithInPoints)
            {
                throw new ApplicationException(string.Format(
                    "the bandwidth must be large enough to accomodate at least {0} points. There are {1} data points, and bandwidth must be at least {2} but it is only {3}",
                    leastBandwithInPoints, n, (double)leastBandwithInPoints / n, bandwidth
                ));
            }

            double[] fitted = null;
            double[] residuals = null;
            double[] robustnessWeights = null;

            if (robustnessIters > 0) { computeFitted = true; } // override computeFitted
            if (computeFitted || robustnessIters > 0) { fitted = new double[n]; }

            // Do an initial fit and 'robustnessIters' robustness iterations.
            // This is equivalent to doing 'robustnessIters+1' robustness iterations
            // starting with all robustness weights set to 1.
            if (robustnessIters > 0)
            {
                residuals = new double[n];
                robustnessWeights = new double[n];
                for (int i = 0; i < robustnessWeights.Length; i++) robustnessWeights[i] = 1;
            }

            for (int iter = 0; iter <= robustnessIters; ++iter)
            {
                int[] bandwidthInterval = { 0, bandwidthInPoints - 1 };
                // At each x, compute a local weighted linear regression
                for (int i = 0; i < n; ++i)
                {
                    double x = xvals[i];

                    // Find out the interval of source points on which
                    // a regression is to be made.
                    if (i > 0)
                    {
                        var newBandwidthInterval = updateBandwidthInterval(x, xvals, bandwidthInterval);
                        if (newBandwidthInterval != null) { bandwidthInterval = newBandwidthInterval; }
                    }

                    // Linear regression to get the coeeficients
                    if (computeFitted || robustnessIters > 0)
                    {
                        double[] coefficients = computeCoefficients(x, xvals, yvals, robustnessWeights, bandwidthInterval);
                        fitted[i] = predict(x, coefficients);
                    }
                    if (robustnessIters > 0) { residuals[i] = Math.Abs(yvals[i] - fitted[i]); }
                }

                // No need to recompute the robustness weights at the last
                // iteration, they won't be needed anymore
                if (iter == robustnessIters) { break; }

                // Recompute the robustness weights.
                double medianResidual = Utilities.Median(residuals); // Find the median residual

                if (medianResidual == 0) { break; }

                for (int i = 0; i < n; ++i)
                {
                    double arg = residuals[i] / (6 * medianResidual);
                    robustnessWeights[i] = (arg >= 1) ? 0 : Math.Pow(1 - arg * arg, 2);
                }
            }

            List<LoessInterval> intervals = computeIntervals(xvals, xStep, bandwidthInPoints);

            return new LoessModel(xvals, yvals, fitted, robustnessWeights, ascendingOrder, intervals);
        }

        private static List<LoessInterval> computeIntervals(double[] xvals, double xStep, int bandwidthInPoints)
        {
            List<LoessInterval> intervals = new List<LoessInterval>();

            int[] bandwidthInterval = new int[] { 0, bandwidthInPoints - 1 };
            double xMin = double.NegativeInfinity;

            for (double x = xvals[0]; x <= xvals[xvals.Length - 1]; x += xStep)
            {
                var newBandwidthInterval = updateBandwidthInterval(x, xvals, bandwidthInterval);
                if (newBandwidthInterval != null)
                {
                    intervals.Add(new LoessInterval(xMin, x, bandwidthInterval));
                    xMin = x;
                    bandwidthInterval = newBandwidthInterval;
                }
            }
            intervals.Add(new LoessInterval(xMin, double.PositiveInfinity, bandwidthInterval));

            return intervals;
        }

        private static double[] computeCoefficients(double x, double[] xval, double[] yval, double[] robustnessWeights, int[] bandwidthInterval)
        {
            int iLeft = bandwidthInterval[0];
            int iRight = bandwidthInterval[1];

            // Compute the point of the bandwidth interval that is
            // farthest from x
            int edge = (x - xval[iLeft] > xval[iRight] - x) ? iLeft : iRight;

            // Compute a least-squares linear fit weighted by
            // the product of robustness weights and the tricube
            // weight function.
            // See http://en.wikipedia.org/wiki/Linear_regression
            // (section "Univariate linear case")
            // and http://en.wikipedia.org/wiki/Weighted_least_squares
            // (section "Weighted least squares")
            double sumWeights = 0;
            double sumX = 0, sumXSquared = 0, sumY = 0, sumXY = 0;
            double denom = Math.Abs(1.0 / (xval[edge] - x));
            for (int k = iLeft; k <= iRight; ++k)
            {
                double xk = xval[k];
                double yk = yval[k];
                double dist = Math.Abs(x - xk);
                double robustnessWeight = robustnessWeights == null ? 1.0 : robustnessWeights[k];
                double w = tricube(dist * denom) * robustnessWeight;
                double xkw = xk * w;
                sumWeights += w;
                sumX += xkw;
                sumXSquared += xk * xkw;
                sumY += yk * w;
                sumXY += yk * xkw;
            }

            double meanX = sumX / sumWeights;
            double meanY = sumY / sumWeights;
            double meanXY = sumXY / sumWeights;
            double meanXSquared = sumXSquared / sumWeights;

            double beta;
            if (meanXSquared == meanX * meanX)
            {
                beta = 0;
            }
            else
            {
                beta = (meanXY - meanX * meanY) / (meanXSquared - meanX * meanX);
            }

            double alpha = meanY - beta * meanX;

            return new double[] { alpha, beta };
        }

        private static double predict(double x, double[] coefficients)
        {
            double y = 0;
            for (int i = 0; i < coefficients.Length; i++)
            {
                y += Math.Pow(x, i) * coefficients[i];
            }

            return y;
        }

        /// <summary>
        /// Given an Index interval into xval that embraces a certain number of points closest to a value less than x,
        /// update the interval so that it embraces the same number of points closest to x.
        /// </summary>
        /// <param name="x">Assumed to be greater than xval[bandwidthInterval[0]] but may be greater than xval[bandwidthInterval[1]]</param>
        /// <param name="xval">Assumed sorted in ascending order</param>
        /// <param name="bandwidthInterval"></param>
        private static int[] updateBandwidthInterval(double x, double[] xval, int[] bandwidthInterval)
        {
            bool updated = false;

            int leftIndex = bandwidthInterval[0];
            int rightIndex = bandwidthInterval[1];

            // x > xval[rightIndex]
            while (rightIndex < xval.Length - 1 && x > xval[rightIndex])
            {
                leftIndex++;
                rightIndex++;
                updated = true;
            }

            // The right edge should be adjusted if the next point to the right
            // is closer to x than the leftmost point of the current interval
            while (rightIndex < xval.Length - 1 && xval[rightIndex + 1] - x < x - xval[leftIndex])
            {
                leftIndex++;
                rightIndex++;
                updated = true;
            }

            if (updated)
            {
                return new int[] { leftIndex, rightIndex };
            }

            return null;
        }

        /**
         * Compute the 
         * <a href="http://en.wikipedia.org/wiki/Local_regression#Weight_function">tricube</a>
         * weight function
         *
         * @param x the argument
         * @return (1-|x|^3)^3
         */
        private static double tricube(double x)
        {
            double tmp = 1 - x * x * x;
            return tmp * tmp * tmp;
        }

        /**
         * Check that all elements of an array are finite real numbers.
         *
         * @param values the values array
         * @param isAbscissae if true, elements are abscissae otherwise they are ordinatae
         * @throws MathException if one of the values is not
         *         a finite real number
         */
        private static void checkAllFiniteReal(double[] values, bool isAbscissae)
        {
            for (int i = 0; i < values.Length; i++)
            {
                double x = values[i];
                if (Double.IsInfinity(x) || Double.IsNaN(x))
                {
                    string pattern = isAbscissae ?
                            "all abscissae must be finite real numbers, but {0}-th is {1}" :
                            "all ordinatae must be finite real numbers, but {0}-th is {1}";
                    throw new ApplicationException(string.Format(pattern, i, x));
                }
            }
        }

        /// <summary>
        /// Check that elements of the abscissae array are in an increasing order.
        /// Throws MathException if the abscissae array is not in an increasing order.
        /// </summary>
        /// <param name="xval">the abscissae array</param>
        private static void checkIncreasing(double[] xval)
        {
            for (int i = 0; i < xval.Length; ++i)
            {
                if (i >= 1 && xval[i - 1] > xval[i])
                {
                    throw new ApplicationException(
                        string.Format("the abscissae array must be sorted in an increasing order, but the {0}-th element is {1} whereas {2}-th is {3}",
                        i - 1, xval[i - 1], i, xval[i]));
                }
            }
        }

        public class LoessModel
        {
            public double[] SortedXs { get; private set; }
            public double[] SortedYs { get; private set; } // sorted in ascending order by X
            public double[] SortedFitted { get; private set; } // sorted in ascending order by X
            public double[] SortedRobustnessWeights { get; private set; } // sorted in ascending order by X
            public int[] OriginalOrder { get; private set; }
            public List<LoessInterval> LoessIntervals { get; private set; }

            public IEnumerable<double> Xs
            {
                get
                {
                    return OriginalOrder.Select(i => SortedXs[i]);
                }
            }

            public IEnumerable<double> Ys
            {
                get
                {
                    return OriginalOrder.Select(i => SortedYs[i]);
                }
            }

            public IEnumerable<double> Fitted
            {
                get
                {
                    if (SortedFitted == null) { return null; }
                    return OriginalOrder.Select(i => SortedFitted[i]);
                }
            }

            public IEnumerable<double> RobustnessWeights
            {
                get
                {
                    if (SortedRobustnessWeights == null) { return null; }
                    return OriginalOrder.Select(i => SortedRobustnessWeights[i]);
                }
            }

            public LoessModel(double[] sortedXs, double[] sortedYs, double[] sortedFitted, double[] sortedRobustnessWeights,
                int[] ascendingOrder, List<LoessInterval> loessIntervals)
            {
                SortedXs = sortedXs;
                SortedYs = sortedYs;
                SortedFitted = sortedFitted;
                SortedRobustnessWeights = sortedRobustnessWeights;
                LoessIntervals = loessIntervals;

                SetOriginalOrder(ascendingOrder);
            }

            private void SetOriginalOrder(int[] ascendingOrder)
            {
                OriginalOrder = new int[ascendingOrder.Length];
                for (int i = 0; i < ascendingOrder.Length; i++)
                {
                    OriginalOrder[ascendingOrder[i]] = i;
                }
            }

            public double[] Predict(IEnumerable<double> xvals)
            {
                double[] xArr = xvals.ToArray();
                int[] ascendingOrder = Enumerable.Range(0, xArr.Length).OrderBy(i => xArr[i]).ToArray(); // argsort xvals
                double[] yArr = new double[xArr.Length];

                int intervalIndex = 0;
                for (int i = 0; i < xArr.Length; i++)
                {
                    double x = xArr[ascendingOrder[i]];
                    while (intervalIndex < LoessIntervals.Count - 1 && LoessIntervals[intervalIndex].IsLeftOf(x))
                    {
                        intervalIndex++;
                    }
                    if (intervalIndex >= LoessIntervals.Count && LoessIntervals[intervalIndex].IsRightOf(x))
                    {
                        throw new ApplicationException(String.Format("Unable to find an interval for {0}.", x));
                    }
                    yArr[ascendingOrder[i]] = predict(x, LoessIntervals[intervalIndex]);
                }

                return yArr;
            }

            public double Predict(double x)
            {
                if (LoessIntervals == null || !LoessIntervals.Any())
                {
                    throw new ApplicationException("Unable to predict because no intervals are available.");
                }

                LoessInterval interval = FindInterval(x);
                if (interval == null)
                {
                    throw new ApplicationException(String.Format("Unable to find an interval for {0}.", x));
                }

                return predict(x, interval);
            }

            private LoessInterval FindInterval(double x)
            {
                int iLeft = 0, iRight = LoessIntervals.Count - 1;
                while (iLeft <= iRight)
                {
                    int iMid = (iLeft + iRight) / 2;
                    if (LoessIntervals[iMid].IsRightOf(x))
                    {
                        iRight = iMid - 1;
                    }
                    else if (LoessIntervals[iMid].IsLeftOf(x))
                    {
                        iLeft = iMid + 1;
                    }
                    else
                    {
                        return LoessIntervals[iMid];
                    }
                }

                return null;
            }

            private double predict(double x, LoessInterval interval)
            {
                var coeffs = LoessInterpolator.computeCoefficients(x, SortedXs, SortedYs, SortedRobustnessWeights, interval.BandwidthInterval);

                return LoessInterpolator.predict(x, coeffs);

            }
        }

        public class LoessInterval
        {
            public double XMin { get; private set; } // inclusive
            public double XMax { get; private set; } // exclusive
            public int[] BandwidthInterval { get; private set; }

            public LoessInterval(double xMin, double xMax, int[] bandwidthInterval)
            {
                XMin = xMin;
                XMax = xMax;
                BandwidthInterval = new int[bandwidthInterval.Length];
                Buffer.BlockCopy(bandwidthInterval, 0, BandwidthInterval, 0, bandwidthInterval.Length * sizeof(int));
            }

            /// <summary>
            /// Whether the interval contains x
            /// </summary>
            /// <param name="x"></param>
            /// <returns></returns>
            public bool Contains(double x)
            {
                return XMin <= x && x < XMax;
            }

            /// <summary>
            /// Whether the interval is on the left of x
            /// </summary>
            /// <param name="x"></param>
            /// <returns></returns>
            public bool IsLeftOf(double x)
            {
                return XMax <= x;
            }

            /// <summary>
            /// Whether the interval is on the right of x
            /// </summary>
            /// <param name="x"></param>
            /// <returns></returns>
            public bool IsRightOf(double x)
            {
                return x < XMin;
            }
        }
    }
}