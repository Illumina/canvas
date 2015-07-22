using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using Illumina.NumericalAnalysis;

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
        private static void GetInnerProdIter(double[] x, double[] I_prod)
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
        private static int GetInnerProdMax(double[] x)
        {
            double[] ipi = new double[x.Length - 1];
            GetInnerProdIter(x, ipi);
            double max_abs_ipi = Math.Abs(ipi[0]);

            for (uint i = 0; i < ipi.Length; i++)
            {
                if (Math.Abs(ipi[i]) > max_abs_ipi)
                {
                    max_abs_ipi = Math.Abs(ipi[i]);
                }
            }
            List<int> indexVector = new List<int>();
            for (int i = 0; i < ipi.Length; i++)
            {
                if (Math.Abs(ipi[i]) == max_abs_ipi)
                {
                    indexVector.Add(i + 1);
                }
            }
            int medIndex = CanvasCommon.Utilities.Median(indexVector);
            return medIndex;
        }

        /// <summary>
        /// given sigma, set Unbalanced Haar wavelet coefficients to zero
        /// </summary>
        private static void HardThresh(List<List<double>> tree, double sigma)
        {
            int J = tree.Count();
            double n = tree[0][5 - 1];
            for (int j = 0; j < J; j++)
            {
                int K = (int)Math.Floor(tree[j].Count() / 5.0);
                for (int k = 0; k < K; k++)
                {
                    if (Math.Abs(tree[j][k * 5 + 2 - 1]) <= sigma * Math.Sqrt(2 * Math.Log(n)))
                    {
                        tree[j][k * 5 + 2 - 1] = 0;
                    }
                }
            } //for (int j=0;j<J;j++)
        }

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
            double[] rec = Enumerable.Repeat(1.0 / Math.Sqrt(n) * smooth, n).ToArray();

            for (int j = 0; j < J; j++)
            {
                int K = (int)Math.Floor((double)tree[j].Count() / 5);
                for (int k = 0; k < K; k++)
                {
                    int skip = k * 5 + 3 - 1;
                    int take = k * 5 + 5 - skip;
                    double[] branch = new double[take];
                    Array.Copy(tree[j].ToArray(), skip, branch, 0, take);
                    List<double> unbal_haar_vector = GetUnbalHaarVector(branch);

                    // for (unsigned i=3; i<=5;i++)
                    for (int i = (int)tree[j][3 - 1 + k * 5] - 1; i < tree[j][5 - 1 + k * 5]; i++)
                    {
                        rec[i] = rec[i] + unbal_haar_vector[i - (int)tree[j][3 - 1 + k * 5] + 1] * tree[j][k * 5 + 2 - 1];
                    }
                }
            }
            return (rec);
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
        /// Find best Unbalanced Haar wavelets decomposition
        /// </summary>
        private static void FindBestUnbalancedHaarDecomposition(double[] x, List<List<double>> tree, double smooth)
        {
            ulong n = (ulong)x.Length;
            tree.Clear();
            List<double> branch = new List<double>();
            branch.Add(0);
            branch.Add(0);
            branch.Add(0);
            branch.Add(0);
            branch.Add(0);

            branch[0] = 1;
            double[] ipi = new double[x.Length - 1];
            GetInnerProdIter(x, ipi);
            int ind_max = GetInnerProdMax(x);
            branch[2] = 1;
            branch[3] = ind_max;
            branch[4] = n;
            branch[1] = ipi[ind_max - 1];
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
                        GetInnerProdIter(subX, ipi);
                        ind_max = GetInnerProdMax(subX);

                        tree[j + 1][(no_child_coeffs - 1) * 5 + 2 - 1] = ipi[ind_max - 1];
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
                        GetInnerProdIter(subX, ipi);
                        ind_max = GetInnerProdMax(subX);

                        tree[j + 1][(no_child_coeffs - 1) * 5 + 2 - 1] = ipi[ind_max - 1];
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

            smooth = 0;
            for (uint i = 0; i < n; i++) smooth += x[i];
            smooth = smooth / Math.Sqrt(n);
        }

        /// <summary>
        /// Enrty function for performing Unbalanced Haar (UH) wavelets decomposition
        /// for change point (breakpoint) detection
        /// </summary>
        public static void HaarWavelets(double[] ratio, double thresholdlower, double thresholdupper, List<int> breakpoints)
        {
            // tree = A list of J matrices (list of lists in C#), where J represents the number of UH “scales”. 
            // Each matrix is ofsize 5 x (the number of UH coefficients at a given scale). 
            // Each column (= vector // of length 5) contains an Unbalanced Haar coefficient in the following format:
            // 1st component - an index of the coefficient; 2nd component - the value of the
            // coefficient; 3rd component - time point where the corresponding UH vector
            // starts; 4th component - last time point before the breakpoint of the UH vector;
            // 5th component - end point of the UH vector.
            List<List<double>> tree = new List<List<double>>();

            double smooth = 0;
            FindBestUnbalancedHaarDecomposition(ratio, tree, smooth);
            // set threshold proportional to SD as and sample size as suggested in
            //  Piotr Fryzlewicz, Journal of the American Statistical Association
            // Vol. 102, No. 480 (Dec., 2007), pp. 1318-1327
            double sd = CanvasCommon.Utilities.StandardDeviation(ratio) * 2.0;
            if (sd < thresholdlower)
            {
                sd = thresholdlower;
            }
            if (sd > thresholdupper)
            {
                sd = thresholdupper;
            }
            HardThresh(tree, sd);
            GetSegments(tree, smooth, breakpoints);


        }
    }
}
