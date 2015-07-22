using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics;

namespace CanvasPartition
{
    /// <summary>
    /// Contains functions/subroutines found in R/getbdry.R and src/getbdry.f
    /// </summary>
    public static class GetBoundary
    {
        /// <summary>
        /// Ported from Fortran
        /// </summary>
        /// <param name="eta"></param>
        /// <param name="nPerm"></param>
        /// <param name="maxOnes"></param>
        /// <param name="sbdry"></param>
        /// <param name="tol"></param>
        public static void ComputeBoundary(double eta, uint nPerm, uint maxOnes, out uint[] sbdry,
            double tol = 1E-2)
        {
            sbdry = new uint[maxOnes * (maxOnes + 1) / 2];
            double[] etaStar = new double[maxOnes];
            double eta0, etaLo, etaHi, pLo, pHi, pExcd;
            uint j, l;

            l = 0;
            sbdry[0] = nPerm - (uint)(nPerm * eta);
            etaStar[0] = eta;
            eta0 = eta;
            for (j = 2; j <= maxOnes; j++)
            {
                etaHi = eta0 * 1.1;
                GetBoundary.EtaBoundary(nPerm, etaHi, j, sbdry, l + 1);
                GetBoundary.PExceed(nPerm, j, sbdry, l + 1, out pHi);
                etaLo = eta0 * 0.25;
                GetBoundary.EtaBoundary(nPerm, etaLo, j, sbdry, l + 1);
                GetBoundary.PExceed(nPerm, j, sbdry, l + 1, out pLo);
                while ((etaHi - etaLo) / etaLo > tol)
                {
                    eta0 = etaLo + (etaHi - etaLo) * (eta - pLo) / (pHi - pLo);
                    GetBoundary.EtaBoundary(nPerm, eta0, j, sbdry, l + 1);
                    GetBoundary.PExceed(nPerm, j, sbdry, l + 1, out pExcd);
                    if (pExcd > eta)
                    {
                        etaHi = eta0;
                        pHi = pExcd;
                    }
                    else
                    {
                        etaLo = eta0;
                        pLo = pExcd;
                    }
                }
                etaStar[j - 1] = eta0;
                l += j;
            }
        }

        /// <summary>
        /// Ported from Fortran
        /// Modifies elements of sbdry
        /// </summary>
        /// <param name="nPerm"></param>
        /// <param name="eta0"></param>
        /// <param name="n1s"></param>
        /// <param name="sbdry"></param>
        /// <param name="sbdryOffset"></param>
        private static void EtaBoundary(uint nPerm, double eta0, uint n1s, uint[] sbdry,
            uint sbdryOffset)
        {
            uint i, k;
            //double dn;
            double tProb;

            double dn = nPerm - n1s; // Hypergeometric parameterizes differently
            k = 0;
            for (i = 1; i <= nPerm; i++)
            {
                // Hypergeometric(total # of balls, # white balls, # balls drawn)
                //var hyper = new Hypergeometric(Convert.ToInt32(nPerm), Convert.ToInt32(n1s),
                //    Convert.ToInt32(i)); // Hypergeometric distribution
                // P(X < k + 0.1) == P(X <= k): k = # of white balls drawn
                //tProb = hyper.CumulativeDistribution(k + 0.1);
                //Console.WriteLine("Hypergeometric({0}, {1}, {2}, {3}) = {4}", nPerm, n1s, i, k, tProb);
                tProb = R.phyper(k, n1s, dn, i, true, false);
                //Console.WriteLine("phyper({0}, {1}, {2}, {3}) = {4}", k, n1s, dn, i, tProb);
                //var hyper = new HypergeometricDistribution(Convert.ToInt32(n1s),
                //    Convert.ToInt32(dn), Convert.ToInt32(i));
                //tProb = hyper.DistributionFunction(Convert.ToInt32(k));
                if (tProb <= eta0)
                {
                    sbdry[sbdryOffset + k] = i;
                    k += 1;
                }
            }
        }

        /// <summary>
        /// Ported from Fortran
        /// Modifies pExcd
        /// </summary>
        /// <param name="nPerm"></param>
        /// <param name="n1s"></param>
        /// <param name="sbdry"></param>
        /// <param name="sbdryOffset"></param>
        /// <param name="pExcd"></param>
        private static void PExceed(uint nPerm, uint n1s, uint[] sbdry, uint sbdryOffset,
            out double pExcd)
        {
            double dlcnk;
            int i, n, k, n1, k1, n2, k2, n3, k3;

            n = Convert.ToInt32(nPerm);
            k = Convert.ToInt32(n1s);
            n1 = Convert.ToInt32(nPerm - sbdry[sbdryOffset]);
            dlcnk = SpecialFunctions.BinomialLn(n, k);
            pExcd = Math.Exp(
                SpecialFunctions.BinomialLn(n1, k) - dlcnk);
            if (n1s >= 2)
            {
                n1 = Convert.ToInt32(sbdry[sbdryOffset]);
                n = Convert.ToInt32(nPerm - sbdry[sbdryOffset + 1]);
                k = Convert.ToInt32(n1s - 1);
                pExcd += Math.Exp(Math.Log(n1) + SpecialFunctions.BinomialLn(n, k) -
                    dlcnk);
            }
            if (n1s >= 3)
            {
                n1 = Convert.ToInt32(sbdry[sbdryOffset]);
                n2 = Convert.ToInt32(sbdry[sbdryOffset + 1]);
                n = Convert.ToInt32(nPerm - sbdry[sbdryOffset + 2]);
                k = Convert.ToInt32(n1s - 2);
                pExcd += Math.Exp(Math.Log(n1) + Math.Log(n1 - 1.0) - Math.Log(2.0) +
                    SpecialFunctions.BinomialLn(n, k) - dlcnk) +
                    Math.Exp(Math.Log(n1) + Math.Log(n2 - n1) +
                    SpecialFunctions.BinomialLn(n, k) - dlcnk);
            }
            if (n1s > 3)
            {
                for (i = 4; i <= n1s; i++)
                {
                    n1 = Convert.ToInt32(sbdry[sbdryOffset + i - 4]);
                    k1 = i - 1;
                    k2 = i - 2;
                    k3 = i - 3;
                    n2 = Convert.ToInt32(sbdry[sbdryOffset + i - 3]);
                    n3 = Convert.ToInt32(sbdry[sbdryOffset + i - 2]);
                    n = Convert.ToInt32(nPerm - sbdry[sbdryOffset + i - 1]);
                    k = Convert.ToInt32(n1s - i + 1);
                    pExcd += Math.Exp(SpecialFunctions.BinomialLn(n1, k1) +
                        SpecialFunctions.BinomialLn(n, k) - dlcnk) +
                        Math.Exp(SpecialFunctions.BinomialLn(n1, k2) +
                        Math.Log(n3 - n1) +
                        SpecialFunctions.BinomialLn(n, k) - dlcnk) +
                        Math.Exp(SpecialFunctions.BinomialLn(n1, k3) +
                        Math.Log(n2 - n1) + Math.Log(n3 - n2) +
                        SpecialFunctions.BinomialLn(n, k) - dlcnk) +
                        Math.Exp(SpecialFunctions.BinomialLn(n1, k3) +
                        Math.Log(n2 - n1) - Math.Log(2.0) + Math.Log(n2 - n1 - 1.0) +
                        SpecialFunctions.BinomialLn(n, k) - dlcnk);

                }
            }
        }

        public static void ComputeBoundary(uint nPerm, double alpha, double eta, out uint[] sbdry)
        {
            GetBoundary.ComputeBoundary(eta, nPerm, Convert.ToUInt32(Math.Floor(nPerm * alpha) + 1), out sbdry);
        }
    }
}
