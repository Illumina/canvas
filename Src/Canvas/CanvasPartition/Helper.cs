using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CanvasPartition
{
    public static class Helper // Helper functions
    {
        /// <summary>
        /// R function: mad.
        /// Assumes that x has no NaN, -Inf, Inf.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="constant"></param>
        /// <returns></returns>
        public static double MedianAbsoluteDeviation(double[] x, double constant = 1.4826)
        {
            double center = Helper.Median(x);
            double[] absDev = new double[x.Length];
            for (int i = 0; i < x.Length; i++) { absDev[i] = Math.Abs(x[i] - center); }
            return constant * Helper.Median(absDev);
        }

        /// <summary>
        /// Computes the median of x[iStart], ..., x[iEnd - 1]
        /// </summary>
        /// <param name="x">the array to compute median over</param>
        /// <param name="iStart">start index: 0-based, inclusive.</param>
        /// <param name="iEnd">end index: 0-based, exclusive; if iEnd == -1, iEnd = x.Length</param>
        /// <returns></returns>
        public static double Median(double[] x, int iStart = 0, int iEnd = -1)
        {
            if (x == null) { return Double.NaN; }
            if (iEnd == -1) { iEnd = x.Length; }

            double[] y = new double[iEnd - iStart]; // Make a copy for QuickSelect
            Array.Copy(x, iStart, y, 0, y.Length);
            int mid = y.Length / 2;
            double median = Helper.QuickSelect(y, mid);
            if (y.Length % 2 == 0) // even
            {
                median = (median + Helper.QuickSelect(y, mid - 1)) / 2;
            }
            return median;
        }

        /// <summary>
        /// Select the k-th smallest item in array x.
        /// Note: it modifies x.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="k">0-based; k == 0 returns the minimum</param>
        /// <returns></returns>
        public static double QuickSelect(double[] x, int k)
        {
            int left = 0, right = x.Length - 1;
            int pos, i;
            double pivot;

            while (left < right)
            {
                pivot = x[k];
                Helper.Swap<double>(ref x[k], ref x[right]);
                for (i = pos = left; i < right; i++)
                {
                    if (x[i] < pivot)
                    {
                        Helper.Swap<double>(ref x[i], ref x[pos]);
                        pos++;
                    }
                }
                Helper.Swap<double>(ref x[right], ref x[pos]);
                if (pos == k) { break; }
                if (pos < k)
                {
                    left = pos + 1;
                }
                else
                {
                    right = pos - 1;
                }
            }
            return x[k];
        }

        /// <summary>
        /// Swap the values of two variables
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="lhs"></param>
        /// <param name="rhs"></param>
        public static void Swap<T>(ref T lhs, ref T rhs)
        {
            T temp;
            temp = lhs;
            lhs = rhs;
            rhs = temp;
        }

        /// <summary>
        /// R function diff: http://stat.ethz.ch/R-manual/R-devel/library/base/html/diff.html
        /// Note: this is a partial implementation
        /// Computes x[i + lag] - x[i]
        /// </summary>
        /// <param name="x"></param>
        /// <param name="lag"></param>
        /// <returns>the difference array</returns>
        public static double[] Diff(double[] x, uint lag = 1)//, uint differences = 1) 
        {
            double[] diff = new double[x.Length - lag];
            for (uint i = lag, j = 0; i < x.Length; i++, j++)
            {
                diff[j] = x[i] - x[j];
            }
            return diff;
        }

        /// <summary>
        /// R function diff: http://stat.ethz.ch/R-manual/R-devel/library/base/html/diff.html
        /// Note: this is a partial implementation
        /// Computes x[i + lag] - x[i]
        /// </summary>
        /// <param name="x"></param>
        /// <param name="lag"></param>
        /// <returns> the difference array</returns>
        public static int[] Diff(int[] x, uint lag = 1)//, uint differences = 1) 
        {
            int[] diff = new int[x.Length - lag];
            for (uint i = lag, j = 0; i < x.Length; i++, j++)
            {
                diff[j] = x[i] - x[j];
            }
            return diff;
        }

        /// <summary>
        /// Computes (sum_i x[i] * w[i]) / (sum_i w[i])
        /// </summary>
        /// <param name="x"></param>
        /// <param name="w"></param>
        /// <param name="iStart">start index: 0-based, inclusive</param>
        /// <param name="iEnd">end index: 0-based, exclusive</param>
        /// <returns></returns>
        public static double WeightedAverage(double[] x, double[] w, int iStart = 0, int iEnd = -1)
        {
            if (iEnd == -1) { iEnd = x.Length; }
            double sumWeight = 0.0;
            double weightedSum = 0.0;
            for (int i = iStart; i < iEnd; i++)
            {
                double wi = (w == null) ? 1.0 : w[i];
                sumWeight += wi;
                weightedSum += x[i] * wi;
            }
            return weightedSum / sumWeight;
        }

        /// <summary>
        /// Computes x[i] = x[i] - y
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static double[] InplaceSub(double[] x, double y)
        {
            for (int i = 0; i < x.Length; i++) { x[i] -= y; }
            return x;
        }

        /// <summary>
        /// Computes x[i] = x[i] - y
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public static int[] InplaceSub(int[] x, int y)
        {
            for (int i = 0; i < x.Length; i++) { x[i] -= y; }
            return x;
        }

        /// <summary>
        /// Computes sum_i w[i] * x[i] * x[i]
        /// </summary>
        /// <param name="x"></param>
        /// <param name="w">w[i] = 1 if w == null</param>
        /// <returns></returns>
        public static double WeightedSumOfSquares(double[] x, double[] w)
        {
            double wss = 0.0;
            for (int i = 0; i < x.Length; i++)
            {
                double wi = (w == null) ? 1.0 : w[i];
                wss += wi * x[i] * x[i];
            }
            return wss;
        }

        /// <summary>
        /// Computes sum_{i=start}^{i=start+size-1} x[i]^pow
        /// </summary>
        /// <param name="x"></param>
        /// <param name="pow"></param>
        /// <param name="start"></param>
        /// <param name="size"></param>
        /// <returns></returns>
        public static double PartialSumOfPowers(double[] x, double pow, int start, int size)
        {
            double sp = 0.0;
            for (int i = start; i < start + size; i++) { sp += Math.Pow(x[i], pow); }
            return sp;
        }

        /// <summary>
        /// Compute cumulative/partial sum.
        /// y[j] = sum_{i=1}^{j} x[i]
        /// </summary>
        /// <param name="x"></param>
        /// <returns>y</returns>
        public static int[] CumulativeSum(int[] x)
        {
            if (x == null || x.Length <= 0) { return null; }
            int[] cumSum = new int[x.Length];
            cumSum[0] = x[0];
            for (int i = 1; i < x.Length; i++) { cumSum[i] = cumSum[i - 1] + x[i]; }
            return cumSum;
        }

        /// <summary>
        /// x[i] = Computes Math.Abs(x[i])
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double[] InplaceAbs(double[] x)
        {
            for (int i = 0; i < x.Length; i++) { x[i] = Math.Abs(x[i]); }
            return x;
        }

        /// <summary>
        /// Computes min_i x[i] and argmin_i x[i] simultaneously
        /// </summary>
        /// <param name="x"></param>
        /// <param name="min">min_i x[i]</param>
        /// <returns>argmin_i x[i]</returns>
        public static int ArgMin(double[] x, out double min)
        {
            if (x == null)
            {
                min = Double.NaN;
                return -1;
            }
            min = x[0];
            int iMin = 0;
            for (int i = 1; i < x.Length; i++)
            {
                if (x[i] < min)
                {
                    min = x[i];
                    iMin = i;
                }
            }
            return iMin;
        }

        /// <summary>
        /// Sort v[i] to v[j-1] in ascending order
        /// Place 1-based argsort indices in indices[i] to indices[j-1]
        /// The FORTRAN interface routines for sorting double precision
        /// vectors are `qsort3' and `qsort4', equivalent to `R_qsort' and
        /// `R_qsort_I', respectively.
        /// — Function: void R_qsort (double *v, int i, int j) 
        /// — Function: void R_qsort_I (double *v, int *I, int i, int j)
        ///  sort v[i:j] (using 1-indexing, i.e., v[1] is the first element)
        /// The ..._I() versions also return the sort.index() vector in I
        /// http://phoenix.inf.upol.cz/cgi-bin/info2www?(R-exts)Utility+functions
        /// </summary>
        /// <param name="v">array to be sorted in-place</param>
        /// <param name="indices">will contain 1-based argsort indices</param>
        /// <param name="i">start index, 0-based inclusive</param>
        /// <param name="j">end index, 0-based exclusive</param>
        public static void QuickSort(double[] v, int[] indices, int i, int j)
        {
            // Uses Array.Sort: not necessarily the quicksort algorithm
            // i: 0-based inclusive
            // j: 0-based exclusive
            int size = j - i;
            if (indices == null)
            {
                Array.Sort<double>(v, i, size);
            }
            else
            {
                for (int k = i; k < j; k++) { indices[k] = k + 1; } // 1-based indices
                Array.Sort<double, int>(v, indices, i, size);
            }
        }

        /// <summary>
        /// R function seq
        /// </summary>
        /// <param name="from"></param>
        /// <param name="to"></param>
        /// <param name="length"></param>
        /// <returns></returns>
        public static double[] Seq(double from, double to, uint length)
        {
            if (length == 0) { return null; }
            double[] seq = new double[length];
            double step = (to - from) / (length - 1);
            seq[0] = from;
            seq[length - 1] = to;
            for (uint i = 1; i < length - 1; i++) { seq[i] = seq[i - 1] + step; }
            return seq;
        }

        /// <summary>
        /// outArray = inArray[indices]
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="inArray"></param>
        /// <param name="indices"></param>
        /// <param name="outArray"></param>
        public static void ExtractValues<T>(T[] inArray, int[] indices, out T[] outArray)
        {
            outArray = new T[indices.Length];
            for (int i = 0; i < indices.Length; i++)
            {
                int j = indices[i];
                outArray[i] = inArray[j];
            }
        }

        /// <summary>
        /// Get indices where the values are not NaN, -Inf, Inf
        /// </summary>
        /// <param name="scores"></param>
        /// <param name="indices"></param>
        public static void GetFiniteIndices(double[] scores, out int[] indices)
        {
            List<int> indexList = new List<int>();
            for (int i = 0; i < scores.Length; i++)
            {
                if (Double.IsInfinity(scores[i]) || Double.IsNaN(scores[i])) { continue; }
                indexList.Add(i);
            }
            indices = indexList.ToArray();
        }
    }
}
