using System;

namespace Illumina.NumericalAnalysis
{
    public static class Selection
    {
        #region members

        private static readonly Random Rand;

        #endregion

        static Selection()
        {
            Rand = new Random();
        }

        /// <summary>
        ///     Uses the quick select algorithm to return the nth element in a list of doubles
        /// </summary>
        public static double QuickSelect(double[] l, int nthElement)
        {
            // sanity check for empty lists
            if (l.Length == 0) return double.NaN;

            return QuickSelectRecursive(l, 0, l.Length - 1, nthElement);
        }

        /// <summary>
        ///     Uses the quick select algorithm to return the nth element in a list of doubles
        /// </summary>
        private static double QuickSelectRecursive(double[] l, int beginIndex, int endIndex, int nthElement)
        {
            // sanity check for lists containing one element
            if (beginIndex == endIndex) return l[beginIndex];

            // pick the pivot
            // N.B. By choosing the pivot point smartly, we can speed up quick select by 1.5x - 2.3x
            int n = endIndex - beginIndex + 1;
            const int theta = 5;

            double pivot = (n <= theta
                                ? l[Rand.Next(beginIndex, endIndex)]
                                : QuickSelectRecursive(l, beginIndex, beginIndex + theta - 1, (theta + 1)/2));

            // partition the data
            int begin = beginIndex;
            int end = endIndex;

            while (begin < end)
            {
                while (l[begin] < pivot) ++begin;
                while (l[end] > pivot) --end;

                if (l[begin] == l[end]) ++begin;
                else if (begin < end)
                {
                    double tmp = l[begin];
                    l[begin] = l[end];
                    l[end] = tmp;
                }
            }

            // shortcut: if begin, end, and n are all the same; all of the values are the same
            if ((begin == end) && (begin == (l.Length - 1))) return pivot;

            int j = end;
            int length = j - beginIndex + 1;

            // recursive calls to quick select
            if (length == nthElement)
            {
                return l[j];
            }
            if (nthElement < length)
            {
                return QuickSelectRecursive(l, beginIndex, j - 1, nthElement);
            }
            return QuickSelectRecursive(l, j + 1, endIndex, nthElement - length);
        }
    }
}