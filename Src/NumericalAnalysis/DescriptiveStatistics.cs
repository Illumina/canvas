using System.Collections.Generic;
using System;

namespace Illumina.NumericalAnalysis
{
    /// <summary>
    ///     This class contains functions that operate on arrays such as:
    ///     max, min, median, mode, std, var
    /// </summary>
    public static class DescriptiveStatistics
    {
        /// <summary>
        ///     Returns the mean of a list of doubles
        /// </summary>
        public static double Mean(List<double> dList)
        {
            double[] dArray = dList.ToArray();
            return Mean(dArray, 0, dArray.Length - 1);
        }

        /// <summary>
        ///     Returns the mean of an array of doubles
        /// </summary>
        public static double Mean(double[] dArray)
        {
            return Mean(dArray, 0, dArray.Length - 1);
        }

        /// <summary>
        ///     Returns the mean of an array of doubles in the range [beginIndex, endIndex]
        /// </summary>
        public static double Mean(double[] dArray, int beginIndex, int endIndex)
        {
            double sum = 0.0;
            for (int i = beginIndex; i <= endIndex; ++i) sum += dArray[i];
            return sum / (endIndex - beginIndex + 1);
        }

        /// <summary>
        ///     Returns the median of an array of doubles
        /// </summary>
        public static double Median(double[] dArray)
        {
            // sanity check: make sure we have some numbers
            int n = dArray.Length;
            if (n == 0) return double.NaN;

            bool isOdd = ((n & 1) != 0);
            int half = (n + 1) / 2;

            // handle the case where we have an odd number of elements
            if (isOdd) return Selection.QuickSelect(dArray, half);

            // handle the case where we have an even number of elements
            return (Selection.QuickSelect(dArray, half) + Selection.QuickSelect(dArray, half + 1)) / 2.0;
        }

        /// <summary>
        ///     Returns the median of a list of doubles
        /// </summary>
        public static double Median(List<double> dList)
        {
            // N.B. converting a list to an array speeds up the median
            // calculation by 74x
            return Median(dList.ToArray());
        }

        /// <summary>
        /// Returns the the median absolute deviation of an array of doubles
        /// </summary>
        /// <param name="constant">mad * 1.4826 estimates standard deviation</param>
        public static double MedianAbsoluteDeviation(double[] dArray, double median = Double.NaN,
            double constant = 1.4826)
        {
            if (Double.IsNaN(median)) { median = Median(dArray); }
            double[] absoluteDeviation = new double[dArray.Length];

            for (int i = 0; i < dArray.Length; i++)
            {
                absoluteDeviation[i] = Math.Abs(dArray[i] - median);
            }
            return Median(absoluteDeviation) * constant;
        }

        /// <summary>
        ///     Returns the the standard deviation of an array of doubles
        /// </summary>
        public static double StandardDeviation(double[] dArray, double mu = Double.NaN)
        {
            double std = 0;

            if (double.IsNaN(mu)) { mu = Mean(dArray); }
            for (int i = 0; i < dArray.Length; i++)
            {
                std += Math.Pow(dArray[i] - mu, 2);
            }
            std = Math.Sqrt(std / dArray.Length);

            return std;
        }

    }
}