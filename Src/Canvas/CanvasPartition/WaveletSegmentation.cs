using System;
using System.Collections.Generic;
using System.Linq;

namespace CanvasPartition
{
    /// <summary>
    /// WavevletSegmentation: Unbalanced HAAR wavelets segmentation 
    /// C# implementation of "Unbalanced Haar Technique for Nonparametric Function Estimation"
    /// Piotr Fryzlewicz, Journal of the American Statistical Association
    /// Vol. 102, No. 480 (Dec., 2007), pp. 1318-1327
    /// </summary>
    internal static class WaveletSegmentation
    {
        /// <summary>
        /// For an input vector of length n, get_inner_prod_iter function computes inner products between the input vector and
        /// all possible n-1 Unbalanced Haar vectors of length n
        /// </summary>
        private static void GetInnerProdIter(double[] x, double[] I_prod, out double mean)
        {
            long n = x.Length;
            double[] I_plus = new double[n - 1];
            double[] I_minus = new double[n - 1];

            I_plus[0] = Math.Sqrt(1 - 1.0 / n) * x[0];
            double sumX = 0; //summing from 2nd element to the end of vector X
            for (uint i = 1; i < n; i++)
            {
                sumX += x[i];
            }
            mean = (x[0] + sumX) / n;
            I_minus[0] = (1.0 / Math.Sqrt(n * (n - 1))) * sumX;

            if (n > 2)
            {
                for (uint m = 1; m < n - 1; m++)
                {
                    double factor = Math.Sqrt((double)((n - m - 1)) * m / (double)(m + 1) / (double)(n - m));
                    I_plus[m] = I_plus[m - 1] * factor + x[m] * Math.Sqrt(1.0 / (m + 1) - 1.0 / n);
                    I_minus[m] = I_minus[m - 1] / factor - (double)(x[m]) / (double)Math.Sqrt(((double)n * n / (double)(m + 1)) - (double)n);
                }
            }
            for (uint i = 0; i < n - 1; i++)
            {
                I_prod[i] = I_plus[i] - I_minus[i];
            }

        }

        /// <summary>
        /// The function finds the Unbalanced Haar vector which yields the largest 
        /// (in absolute value) inner product with the input vector
        /// </summary>
        private static int GetInnerProdMax(double[] ipi)
        {
            double max_abs_ipi = ipi.Select(ip => Math.Abs(ip)).Max();
            int index = 0;
            for (; index < ipi.Length; index++)
            {
                if (Math.Abs(ipi[index]) == max_abs_ipi)
                    break;
            }
            // this function used to return the median of best index values in case of a tie;
            // that seems silly, especially since if there are an even number, the median will 
            // return something between the best solutions!
            return index + 1;
        }

        /// <summary>
        /// given sigma, set Unbalanced Haar wavelet coefficients to zero
        /// </summary>
        private static void HardThresh(List<List<double>> tree, double sigma, bool isGermline, string chromosome)
        {
            int treeSize = tree.Count();
            List<double> thresholds = new List<double>();
            int[] indices = new int[treeSize];
            // makes threshold depend on the complexity of wavelet basis functions 
            if (isGermline)
            {
                int[] counts = new int[treeSize];
                for (int nodeIndex = 0; nodeIndex < treeSize; nodeIndex++)
                {
                    counts[nodeIndex] = (int)Math.Floor(tree[nodeIndex].Count() / 5.0);
                    indices[nodeIndex] = nodeIndex;
                }
                // transform number of segments at each wavelet scale into NewMin-NewMax range, use it to weight threshold function
                Array.Sort<int>(indices, (a, b) => counts[b].CompareTo(counts[a]));
                double NewMax = 1.0;
                double NewMin = 0.8;
                thresholds = Enumerable.Range(1, treeSize).Select(x => ((double)x * (NewMax - NewMin)) / treeSize + NewMin).ToList();
            }
            else
            {
                for (int nodeIndex = 0; nodeIndex < treeSize; nodeIndex++)
                {
                    thresholds.Add(1.0);
                    indices[nodeIndex] = nodeIndex;
                }
            }
            int subtreeSize = 5;
            double n = tree[0][subtreeSize - 1];
            for (int nodeIndex = 0; nodeIndex < treeSize; nodeIndex++)
            {
                int K = (int)Math.Floor(tree[nodeIndex].Count() / 5.0);
                for (int k = 0; k < K; k++)
                {
                    // threshold 2x*sigma for TN and ranges from 0.5x to 1.5x *sigma for Germline
                    // parameters trained on CanvasRegression datasets                     
                    if (Math.Abs(tree[nodeIndex][k * subtreeSize + 2 - 1]) <= 2 * sigma * (thresholds[indices[nodeIndex]]) * Math.Sqrt(2 * Math.Log(n)))
                    {
                        tree[nodeIndex][k * subtreeSize + 2 - 1] = 0;
                    }
                }
            }
        }

        // 
        private static List<double> GetUnbalHaarVector(double[] a)
        {
            double n = a[2] - a[0] + 1;
            double m = a[1] - a[0] + 1;
            List<double> returnVector = new List<double>();
            double val1 = Math.Sqrt(1 / m - 1 / n);
            double val2 = -1.0 / Math.Sqrt(n * n / m - n);

            for (uint i = 0; i < Math.Floor(n); i++)
            {
                if (i < m) returnVector.Add(val1);
                else returnVector.Add(val2);
            }
            return returnVector;
        }

        /// <summary>
        /// Reconstructs a vector from its top-down Unbalanced Haar decomposition stored in an object returned
        /// by best.unbal.haar or hard.thresh
        /// </summary>
        private static double[] GetReconstructedVector(List<List<double>> tree, double smooth)
        {
            int J = tree.Count();

            int n = (int)tree[0][5 - 1];
            double[] rec = new double[n];
            for (int i = 0; i < n; i++)
            {
                rec[i] = 1.0 / Math.Sqrt(n) * smooth;
            }

            for (int j = 0; j < J; j++)
            {
                int K = (int)Math.Floor((double)tree[j].Count() / 5);
                for (int k = 0; k < K; k++)
                {
                    int skip = k * 5 + 3 - 1;
                    int take = k * 5 + 5 - skip;
                    double[] branch = new double[take];
                    for (int i = 0; i < take; i++)
                        branch[i] = tree[j][skip + i];
                    List<double> unbal_haar_vector = GetUnbalHaarVector(branch);

                    for (int i = (int)tree[j][3 - 1 + k * 5] - 1; i < tree[j][5 - 1 + k * 5]; i++)
                    {
                        rec[i] = rec[i] + unbal_haar_vector[i - (int)tree[j][3 - 1 + k * 5] + 1] * tree[j][k * 5 + 2 - 1];
                    }
                }
            }
            return rec;
        }

        /// <summary>
        /// Function outputs breakpoints from a reconstructed vector given 
        /// its top-down Unbalanced Haar decomposition 
        /// </summary>
        private static void GetSegments(List<List<double>> tree, double smooth, List<int> breakpoints)
        {
            double[] recontr_seq = GetReconstructedVector(tree, smooth);
            breakpoints.Add(0);
            for (int i = 1; i < (int)recontr_seq.Length; i++)
            {
                if (recontr_seq[i] - recontr_seq[i - 1] != 0)
                {
                    breakpoints.Add(i);
                }
            }
        }

        /// <summary>
        /// Eliminate unneeded segment boundaries, mostly coarse-grained splits due to waviness;
        /// </summary>
        /// <param name="prelimBreakpoints"></param>
        /// <param name="tree"></param>
        /// <param name="v"></param>
        /// <returns>revised set of breakpoints with poorly-supported breakpoints removed</returns>
        private static void GetBreakpointsAfterHealingBadSplits(List<int> breakpoints, List<int> prelimBreakpoints, double[] ratio, List<double> factorOfThreeCMADs, string chromosome)
        {
            // Wavelet segmentation can introduce breakpoints that are not well-supported. 
            // This happens when something that should be segmented as 0000000001111111112222222222
            // goes through a series of steps that chop it more like this:
            //                 split 1                                 xxxxxxxxxxxxxxxxxxyyyyyyyyyy
            //                 split 2                                 wwwwwxxxxxxxxxxxxxyyyyyyyyyy
            //                 split 3                                 wwwwwzzzzxxxxxxxxxyyyyyyyyyy
            //                                                         
            // If split 2 scores well enough, fundamentally because 11111 != 00000, we will
            // end up including the breakpoint between www and zzz even though the difference between
            // them is small.  This step tests whether the medians of the implied segments are really
            // different enough to retain the breakpoint.
            //
            // The process is greedy, proceeding from left to right.  It might be better to merge
            // most-similar adjacent segments first, but this seems good enough for now.
            int N = ratio.Length;
            int L = prelimBreakpoints.Count();
            breakpoints.Add(prelimBreakpoints[0]);
            for (int i = 1; i < L; ++i)
            {
                var leftStart = breakpoints[breakpoints.Count() - 1];
                var rightStart = prelimBreakpoints[i];
                var rightEnd = (i < L - 1) ? prelimBreakpoints[i + 1] : N;
                var leftLength = rightStart - leftStart;
                var rightLength = rightEnd - rightStart;
                var leftMedian = CanvasCommon.Utilities.Median(ratio.Skip(leftStart).Take(leftLength));
                var rightMedian = CanvasCommon.Utilities.Median(ratio.Skip(rightStart).Take(rightLength));
                double weightedMedian = (leftLength * leftMedian + rightLength * rightMedian) / (rightEnd - leftStart);
                int smallerLength = Math.Min(leftLength,rightLength);
                var factorOfThreeScale = Math.Min(factorOfThreeCMADs.Count() - 1, (int)Math.Ceiling(Math.Log(smallerLength) / Math.Log(3)));
                double factorOfThreeCutoff = factorOfThreeCMADs[factorOfThreeScale];
                //Console.WriteLine($"Healing {chromosome}: {leftStart} {rightStart} {rightEnd} {leftMedian} {rightMedian} {factorOfThreeCutoff}");
                if (Math.Abs(leftMedian - rightMedian) > factorOfThreeCutoff * 4 * Math.Max(weightedMedian, 50d))
                {
                    breakpoints.Add(prelimBreakpoints[i]);
                }
            }
        }

        /// <summary>
        /// Refines breakpoint by maximising local medians
        /// </summary>
        private static void RefineSegments(List<int> breakpoints, IEnumerable<double> coverage)
        {
            int halfWindow = 5;
            double totalMedian = CanvasCommon.Utilities.Median(coverage);
            for (int i = 1; i < breakpoints.Count - 1; i++)
            {
                int leftInterval = Math.Min(halfWindow, (breakpoints[i] - breakpoints[i - 1]) / 2);
                int rightInterval = Math.Min(halfWindow, (breakpoints[i + 1] - breakpoints[i]) / 2);
                double bestMedianDifference = Math.Abs(CanvasCommon.Utilities.Median(coverage, breakpoints[i - 1], breakpoints[i]) - totalMedian);
                int bestBreakpoint = breakpoints[i];
                for (int j = breakpoints[i] - leftInterval; j < breakpoints[i] + rightInterval; j++)
                {
                    double tmpMedianDifference = Math.Abs(CanvasCommon.Utilities.Median(coverage, breakpoints[i - 1], j) - totalMedian);
                    if (tmpMedianDifference > bestMedianDifference)
                    {
                        bestMedianDifference = tmpMedianDifference;
                        bestBreakpoint = j;
                    }
                }
                breakpoints[i] = bestBreakpoint;
            }
        }

        /// <summary>
        /// Find best Unbalanced Haar wavelets decomposition
        /// </summary>
        /// <returns>measure of smoothness (not used?)</returns>
        private static double FindBestUnbalancedHaarDecomposition(double[] x, List<List<double>> tree)
        {
            ulong n = (ulong)x.Length;
            tree.Clear();
            var branch = new List<double>();
            branch.Add(0);
            branch.Add(0);
            branch.Add(0);
            branch.Add(0);
            branch.Add(0);

            branch[0] = 1;
            double[] ipi = new double[x.Length - 1];
            GetInnerProdIter(x, ipi, out double mean);
            int ind_max = GetInnerProdMax(ipi);
            branch[2] = 1;
            branch[3] = ind_max;
            branch[4] = n;
            double meanscale = 200.0;
            branch[1] = ipi[ind_max - 1] / Math.Max(0.5, mean / meanscale);
            tree.Add(new List<double>());
            for (int i = 0; i < (int)branch.Count(); i++)
            {
                tree[0].Add(branch[i]);
            }

            int j = 0;
            double bpSum = 0;

            for (int i = 0; i < (int)Math.Floor((double)tree[j].Count() / 5); i++)
            {
                bpSum += tree[j][5 + i * 5 - 1] - tree[j][3 + i * 5 - 1] - 1.0;
            }

            while (bpSum != 0)
            {
                int no_parent_coeffs = (int)Math.Floor((double)tree[j].Count() / 5);
                int no_child_coeffs = 0;
                for (int i = 0; i < no_parent_coeffs; i++)
                {

                    if (tree[j][4 + 5 * i - 1] - tree[j][3 + 5 * i - 1] >= 1)
                    {
                        no_child_coeffs++;

                        if (tree.Count() == j + 1)
                        {
                            tree.Add(new List<double>() { 0, 0, 0, 0, 0 });
                        }
                        else
                        {
                            for (uint k = 0; k < 5; k++) tree[j + 1].Add(0);
                        }

                        tree[j + 1][1 + 5 * (no_child_coeffs - 1) - 1] = 2 * tree[j][1 + 5 * i - 1] - 1;
                        int skip = (int)tree[j][i * 5 + 3 - 1] - 1;
                        int take = (int)tree[j][i * 5 + 4 - 1] - skip;
                        double[] subX = new double[take];
                        Array.Copy(x, skip, subX, 0, take);
                        Array.Clear(ipi, 0, ipi.Length);
                        Array.Resize(ref ipi, subX.Length - 1);
                        GetInnerProdIter(subX, ipi, out double leftmean);
                        ind_max = GetInnerProdMax(ipi);

                        tree[j + 1][(no_child_coeffs - 1) * 5 + 2 - 1] = ipi[ind_max - 1] / Math.Max(0.5, leftmean / meanscale);
                        tree[j + 1][(no_child_coeffs - 1) * 5 + 3 - 1] = tree[j][i * 5 + 3 - 1];
                        tree[j + 1][(no_child_coeffs - 1) * 5 + 5 - 1] = tree[j][i * 5 + 4 - 1];
                        tree[j + 1][(no_child_coeffs - 1) * 5 + 4 - 1] = ind_max + tree[j][i * 5 + 3 - 1] - 1;

                    } //end if (tree.at(j)[4+5*(i-1)-1]-tree.at(j)[3+5*(i-1)-1]>=1) {


                    if (tree[j][i * 5 + 5 - 1] - tree[j][i * 5 + 4 - 1] >= 2)
                    {
                        no_child_coeffs++;

                        if (tree.Count() == j + 1)
                        {
                            tree.Add(new List<double>() { 0, 0, 0, 0, 0 });
                        }
                        else
                        {
                            for (uint k = 0; k < 5; k++) tree[j + 1].Add(0);
                        }


                        tree[j + 1][1 + 5 * (no_child_coeffs - 1) - 1] = 2 * tree[j][1 + 5 * i - 1];
                        int skip = (int)tree[j][i * 5 + 4 - 1];
                        int take = (int)tree[j][i * 5 + 5 - 1] - skip;
                        double[] subX = new double[take];
                        Array.Copy(x, skip, subX, 0, take);
                        Array.Clear(ipi, 0, ipi.Length);
                        Array.Resize(ref ipi, subX.Length - 1);
                        GetInnerProdIter(subX, ipi, out double rightmean);
                        ind_max = GetInnerProdMax(ipi);

                        tree[j + 1][(no_child_coeffs - 1) * 5 + 2 - 1] = ipi[ind_max - 1] / Math.Max(0.5, rightmean / meanscale);
                        tree[j + 1][(no_child_coeffs - 1) * 5 + 3 - 1] = tree[j][i * 5 + 4 - 1] + 1;
                        tree[j + 1][(no_child_coeffs - 1) * 5 + 5 - 1] = tree[j][i * 5 + 5 - 1];
                        tree[j + 1][(no_child_coeffs - 1) * 5 + 4 - 1] = ind_max + tree[j][i * 5 + 4 - 1];
                    }
                }
                j++;

                bpSum = 0;
                for (int k = 0; k < (int)Math.Floor((double)tree[j].Count() / 5); k++)
                {
                    bpSum += tree[j][5 + k * 5 - 1] - tree[j][3 + k * 5 - 1] - 1;
                }

            } //while 1

            double smooth = 0;
            for (uint i = 0; i < n; i++) smooth += x[i];
            return smooth / Math.Sqrt(n);
        }

        /// <summary>
        /// Entry function for performing Unbalanced Haar (UH) wavelets decomposition
        /// for change point (breakpoint) detection
        /// </summary>
        public static void HaarWavelets(double[] ratio, double thresholdlower, double thresholdupper, List<int> breakpoints,
            bool isGermline, double madFactor, double? coeffVariability, List<double> factorOfThreeCMADs, string chromosome = "")
        {
            // tree = A list of J matrices (list of lists in C#), where J represents the number of UH “scales”. 
            // Each matrix is of size 5 x (the number of UH coefficients at a given scale). 
            // Each column (= vector // of length 5) contains an Unbalanced Haar coefficient in the following format:
            // 1st component - an index of the coefficient; 2nd component - the value of the
            // coefficient; 3rd component - time point where the corresponding UH vector
            // starts; 4th component - last time point before the breakpoint of the UH vector;
            // 5th component - end point of the UH vector.

            var tree = new List<List<double>>();
            double smooth = FindBestUnbalancedHaarDecomposition(ratio, tree);

            // Originally:
            // set threshold proportional to SD as and sample size as suggested in
            // Piotr Fryzlewicz, Journal of the American Statistical Association
            // Vol. 102, No. 480 (Dec., 2007), pp. 1318-1327
            //
            // Now:
            // Set threshold as the coefficient of variation times the median coverage, times the 'madFactor'
            double median = CanvasCommon.Utilities.Median(ratio);
            double variabilityMeasure = coeffVariability.HasValue ? median * coeffVariability.Value : CanvasCommon.Utilities.Mad(ratio);
            double threshold = madFactor * variabilityMeasure;

            if (threshold < thresholdlower)
            {
                threshold = thresholdlower;
            }
            if (threshold > thresholdupper)
            {
                threshold = thresholdupper;
            }

            HardThresh(tree, threshold, isGermline, chromosome);

            var prelimBreakpoints = new List<int> { };
            GetSegments(tree, smooth, prelimBreakpoints);
            GetBreakpointsAfterHealingBadSplits(breakpoints, prelimBreakpoints, ratio, factorOfThreeCMADs, chromosome);
            if (isGermline)
                RefineSegments(breakpoints, ratio);
        }
    }
}
