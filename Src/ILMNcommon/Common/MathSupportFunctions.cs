using System;
using System.Collections.Generic;
using System.Linq;

namespace Illumina.Common
{
    public static class MathSupportFunctions
    {
        public static double Correlate(List<double> listA, List<double> listB)
        {
            if (listA.Count == 0 || listA.Count != listB.Count)
            {
                return 0;
            }
            double meanA = 0;
            double meanB = 0;
            for (int index = 0; index < listA.Count; index++)
            {
                meanA += listA[index];
                meanB += listB[index];
            }
            meanA /= listA.Count;
            meanB /= listB.Count;
            double stdDevA = 0;
            double stdDevB = 0;
            double correlation = 0;
            for (int index = 0; index < listA.Count; index++)
            {
                double diffA = listA[index] - meanA;
                double diffB = listB[index] - meanB;
                stdDevA += (diffA * diffA);
                stdDevB += (diffB * diffB);
                correlation += diffA * diffB;
            }
            stdDevA = (float)Math.Sqrt(stdDevA / listA.Count);
            stdDevB = (float)Math.Sqrt(stdDevB / listA.Count);
            correlation /= (stdDevA * stdDevB * listA.Count);
            return correlation;
        }

        /// <summary>
        ///     Helper function: Compute the best least-squares fit of a line to a set of x,y values.
        ///     Returns a float array containing (offset, slope)
        /// </summary>
        public static float[] LinearFit(float[] x, float[] y)
        {
            float sx = 0; // Sum(x)
            float sy = 0; // Sum(y)
            float sxx = 0; // Sum(x*x)
            float sxy = 0; // Sum(x*y)

            float slope = 0;
            float offset = 0;

            float[] returnValues = new float[]
            {
                0, 0
            };

            if (x != null && y != null)
            {
                if (x.Length == y.Length)
                {
                    for (int i = 0; i < x.Length; i++) sx += x[i];
                    for (int i = 0; i < x.Length; i++) sy += y[i];
                    for (int i = 0; i < x.Length; i++) sxy += x[i] * y[i];
                    for (int i = 0; i < x.Length; i++) sxx += x[i] * x[i];
                    float denominator = (x.Length) * sxx - sx * sx;
                    if (denominator != 0)
                    {
                        slope = (x.Length * sxy - sx * sy) / denominator;
                        offset = (sy * sxx - sx * sxy) / denominator;
                    }
                }
            }
            returnValues[0] = offset;
            returnValues[1] = slope;
            return returnValues;
        }

        public static bool SolveLinearEquation(out double pK1, out double pK2, out double pK3,
            double a, double b, double c, double d, double e, double f,
            double v1, double v2, double v3)
        {
            pK1 = pK2 = pK3 = 0;

            double k1 = d * f - e * e;
            double k2 = c * e - b * f;
            double k3 = b * e - c * d;
            double k4 = a * f - c * c;
            double k5 = c * b - a * e;
            double k6 = a * d - b * b;

            double dd = a * k1 + b * k2 + c * k3;

            if (Math.Abs(dd) < 1e-4)
                return false;

            pK1 = k1 * v1 + k2 * v2 + k3 * v3;
            pK2 = k2 * v1 + k4 * v2 + k5 * v3;
            pK3 = k3 * v1 + k5 * v2 + k6 * v3;

            pK1 /= dd;
            pK2 /= dd;
            pK3 /= dd;

            return true;
        }

        public static float Mean(float[] vec)
        {
            if (vec.Length == 0) return float.NaN;
            double sum = 0;
            foreach (float f in vec) sum += f;
            return (float)(sum / vec.Length);
        }

        public static float Stdev(float[] vec, float meanval)
        {
            if (vec.Length <= 1) return 0;
            double sum = 0;
            foreach (float f in vec)
            {
                float val = f - meanval;
                sum += val * val;
            }
            return (float)Math.Sqrt(sum / (vec.Length - 1));
        }

        public static float Stdev(float[] vec)
        {
            float meanval = Mean(vec);
            return Stdev(vec, meanval);
        }

        public static float Min(float[] vec)
        {
            float minval = float.MaxValue;
            foreach (float f in vec)
            {
                if (float.IsNaN(f)) continue;
                if (f < minval)
                    minval = f;
            }
            return minval;
        }

        public static float Max(float[] vec)
        {
            float maxval = float.MinValue;
            foreach (float f in vec)
            {
                if (float.IsNaN(f)) continue;
                if (f > maxval)
                    maxval = f;
            }
            return maxval;
        }

        public static int Percentile2Index(decimal percentile, int startIndex, int endIndex)
        {
            if (percentile < 0.0m || percentile > 1.0m)
            {
                throw new ArgumentException(String.Format("Asking for the {0}th percentile.", percentile * 100));
            }
            if (startIndex < 0)
            {
                throw new ArgumentException(String.Format("The start index ({0}) is negative.", startIndex));
            }
            if (startIndex > endIndex)
            {
                throw new ArgumentException(String.Format("The startIndex ({0}) must be <= the endIndex ({1})", startIndex, endIndex));
            }

            int numElements = endIndex - startIndex + 1;
            return startIndex + Math.Max(0, ((int)Math.Ceiling(percentile * numElements)) - 1);
        }

        public static float GetPercentileHandleNaN(float[] values, decimal Percentile)
        {
            if (values.Length == 0)
                return float.NaN;

            int nanIndex = 0;
            int maxIndex = values.Length - 1;
            int realIndex = maxIndex;
            int lastRealIndex;

            //swap all the NaNs over to the right
            while (true)
            {
                //walk in from the left and right, to find the first NaN and non-NaN, respectively
                float value = values[realIndex];
#pragma warning disable 1718
                while (value != value)
                {
                    if (realIndex == 0) return float.NaN;
                    value = values[--realIndex];
                }

                value = values[nanIndex];
                while (value == value)
                {
                    if (nanIndex == maxIndex) break;
                    value = values[++nanIndex];
                }
#pragma warning restore 1718
                lastRealIndex = realIndex;

                //make a swap
                if (nanIndex < realIndex)
                {
                    float temp = values[realIndex];
                    values[realIndex] = values[nanIndex];
                    values[nanIndex] = temp;
                }
                else break;
            }

            int PercentileIndex = Percentile2Index(Percentile, 0, lastRealIndex);
            return GetPercentileNoNaNs(values, 0, lastRealIndex, PercentileIndex);
        }

        public static float GetPercentileNoNaNs(float[] values, int rangeStart, int rangeEnd, decimal percentile)
        {
            int percentileIndex = Percentile2Index(percentile, rangeStart, rangeEnd);
            return GetPercentileNoNaNs(values, rangeStart, rangeEnd, percentileIndex);
        }

        public static ushort GetPercentileNoNaNs(ushort[] values, int rangeStart, int rangeEnd, decimal percentile)
        {
            int percentileIndex = Percentile2Index(percentile, rangeStart, rangeEnd);
            return GetPercentileNoNaNs(values, rangeStart, rangeEnd, percentileIndex);
        }

        /// <summary>
        /// Quickly find the TargetIndex-th element in Values[RangeStart...RangeEnd].  
        /// Assumes 0 <= TargetIndex < RangeSize
        /// Adapted from the RANDOMIZED-SELECT pseudocode in "Introduction to Algorithms" by Cormen et al, page 186
        /// (chapter 9.2, Selection in expected linear time)
        /// </summary>
        private static float GetPercentileNoNaNs(float[] values, int rangeStart, int rangeEnd, int targetIndex)
        {
            if (rangeEnd >= values.Length)
                throw new ArgumentException(String.Format("RangeEnd ({0}) is >= Values.Length ({1})", rangeEnd, values.Length));

            Random Randomizer = new Random();

            while (true)
            {
                if (rangeEnd == rangeStart) return values[rangeStart];

                // Randomized partition of Values[RangeStart...RangeEnd]:
                int swapIndex = rangeStart + Randomizer.Next() % (rangeEnd - rangeStart);
                float swapValue = values[swapIndex];
                values[swapIndex] = values[rangeEnd];
                values[rangeEnd] = swapValue;

                int i = rangeStart - 1;
                float pivot = values[rangeEnd];
                int countBelowPivot = 0;
                for (int j = rangeStart; j < rangeEnd; j++)
                {
                    if (values[j] <= pivot)
                    {
                        if (values[j] < pivot) countBelowPivot++;
                        i++;
                        swapValue = values[j];
                        values[j] = values[i];
                        values[i] = swapValue;
                    }
                }

                swapValue = values[i + 1];
                values[i + 1] = values[rangeEnd];
                values[rangeEnd] = swapValue;
                int rangeMiddle = i + 1;
                // Now values from RangeStart...RangeMiddle are <= Pivot, and values from RangeMiddle+1...RangeEnd are > Pivot

                // performance fix for equal values array:
                if (countBelowPivot == 0 && rangeMiddle == rangeEnd)
                    return values[rangeStart];

                if (targetIndex < rangeMiddle)
                {
                    rangeEnd = rangeMiddle - 1;
                }
                else if (targetIndex > rangeMiddle)
                {
                    rangeStart = rangeMiddle + 1;
                }
                else
                {
                    return values[rangeMiddle];
                }
            }
        }

        /// <summary>
        /// Quickly find the TargetIndex-th element in Values[RangeStart...RangeEnd].  
        /// Assumes 0 <= TargetIndex < RangeSize
        /// Adapted from the RANDOMIZED-SELECT pseudocode in "Introduction to Algorithms" by Cormen et al, page 186
        /// (chapter 9.2, Selection in expected linear time)
        /// </summary>
        private static ushort GetPercentileNoNaNs(ushort[] values, int rangeStart, int rangeEnd, int targetIndex)
        {
            if (rangeEnd >= values.Length)
                throw new ArgumentException(String.Format("RangeEnd ({0}) is >= Values.Length ({1})", rangeEnd, values.Length));

            Random Randomizer = new Random();

            while (true)
            {
                if (rangeEnd == rangeStart) return values[rangeStart];

                // Randomized partition of Values[RangeStart...RangeEnd]:
                int swapIndex = rangeStart + Randomizer.Next() % (rangeEnd - rangeStart);
                ushort swapValue = values[swapIndex];
                values[swapIndex] = values[rangeEnd];
                values[rangeEnd] = swapValue;

                int i = rangeStart - 1;
                ushort pivot = values[rangeEnd];
                int countBelowPivot = 0;
                for (int j = rangeStart; j < rangeEnd; j++)
                {
                    if (values[j] <= pivot)
                    {
                        if (values[j] < pivot) countBelowPivot++;
                        i++;
                        swapValue = values[j];
                        values[j] = values[i];
                        values[i] = swapValue;
                    }
                }

                swapValue = values[i + 1];
                values[i + 1] = values[rangeEnd];
                values[rangeEnd] = swapValue;
                int rangeMiddle = i + 1;
                // Now values from RangeStart...RangeMiddle are <= Pivot, and values from RangeMiddle+1...RangeEnd are > Pivot

                // performance fix for equal values array:
                if (countBelowPivot == 0 && rangeMiddle == rangeEnd)
                    return values[rangeStart];

                if (targetIndex < rangeMiddle)
                {
                    rangeEnd = rangeMiddle - 1;
                }
                else if (targetIndex > rangeMiddle)
                {
                    rangeStart = rangeMiddle + 1;
                }
                else
                {
                    return values[rangeMiddle];
                }
            }
        }

        // assumes vals has already been sorted
        // uses linear interpolation
        public static float PercentilePreSorted(float[] vals, int percentile)
        {
            int n = vals.Length;
            if (n == 0)
                return float.NaN;
            int i1 = n * percentile / 100;
            float f = n * percentile / 100.0f - i1;
            if (f < 0.5f) i1--;
            if (i1 < 0)
                return vals[0];
            if (i1 >= n - 1)
                return vals[n - 1];
            float x1 = 100 * (i1 + 0.5f) / n;
            float x2 = 100 * (i1 + 1 + 0.5f) / n;
            float y1 = vals[i1];
            float y2 = vals[i1 + 1];
            float m = (y2 - y1) / (x2 - x1);
            return y1 + m * (percentile - x1);
        }

        public static double Median<TKey>(IEnumerable<KeyValuePair<TKey, long>> sortedHistogram) where TKey : struct
        {
            long totalValues = sortedHistogram.Sum(kvp => kvp.Value);
            long cumulativeValuesCount = 0;
            long lowerMedianValuesCount = (totalValues + 1) / 2;
            bool evenNumberOfValues = totalValues % 2 == 0;
            long upperMedianValuesCount = evenNumberOfValues ? lowerMedianValuesCount + 1 : lowerMedianValuesCount;
            TKey? lowerMedian = null;
            foreach (var histogramEntry in sortedHistogram)
            {
                var value = histogramEntry.Key;
                cumulativeValuesCount += histogramEntry.Value;
                if (!lowerMedian.HasValue && cumulativeValuesCount >= lowerMedianValuesCount)
                    lowerMedian = value;
                if (cumulativeValuesCount >= upperMedianValuesCount)
                    return ((dynamic)lowerMedian + value) / 2.0;
            }
            throw new ApplicationException("Median is undefined");
        }
    }
}