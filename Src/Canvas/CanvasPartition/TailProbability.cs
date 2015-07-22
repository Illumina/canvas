using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.Distributions;

namespace CanvasPartition
{
    /// <summary>
    /// Contains subroutines/functions found in src/tailprobs.f
    /// </summary>
    internal static class TailProbability
    {
        /// <summary>
        /// tail probability of circular binary segmentation statistic
        /// from Siegmund (1988) or Yao (1989) paper
        /// </summary>
        /// <param name="b"></param>
        /// <param name="delta"></param>
        /// <param name="m"></param>
        /// <param name="nGrid"></param>
        /// <param name="tol"></param>
        /// <returns></returns>
        public static double TailP(double b, double delta, int m, int nGrid, double tol)
        {
            double t, tl, dincr, bsqrtm, x, nux;

            dincr = (0.5 - delta) / nGrid;
            bsqrtm = b / Math.Sqrt(m);

            tl = 0.5 - dincr;
            t = 0.5 - 0.5 * dincr;
            double tailP = 0.0;
            for (int i = 0; i < nGrid; i++)
            {
                tl = tl + dincr;
                t = t + dincr;
                x = bsqrtm / Math.Sqrt(t * (1 - t));
                nux = TailProbability.Nu(x, tol);
                tailP = tailP + Math.Pow(nux, 2) * TailProbability.IntegralInvT1tSq(tl, dincr);
            }
            tailP = 9.973557E-2 * Math.Pow(b, 3) * Math.Exp(-Math.Pow(b, 2) / 2) * tailP;
            // since test is two-sided need to multiply tailp by 2
            tailP = 2.0 * tailP;

            return tailP;
        }

        private static double Nu(double x, double tol)
        {
            double nu;
            double lnu0, lnu1, dk, xk;
            int i, k;

            // fpnorm(x): pnorm(x, 0, 1, lower.tail=TRUE, log.p=FALSE)
            // calculates P(X <= x)
            var norm = new Normal(); // N(0, 1)
            if (x > 0.01)
            {
                lnu1 = Math.Log(2.0) - 2 * Math.Log(x);
                lnu0 = lnu1;
                k = 2;
                dk = 0;
                for (i = 0; i < k; i++)
                {
                    dk = dk + 1;
                    xk = -x * Math.Sqrt(dk) / 2.0;
                    lnu1 = lnu1 - 2.0 * norm.CumulativeDistribution(xk) / dk;
                }
                while (Math.Abs((lnu1 - lnu0) / lnu1) > tol)
                {
                    lnu0 = lnu1;
                    for (i = 0; i < k; i++)
                    {
                        dk = dk + 1;
                        xk = -x * Math.Sqrt(dk) / 2.0;
                        lnu1 = lnu1 - 2.0 * norm.CumulativeDistribution(xk) / dk;
                    }
                    k *= 2;
                }
            }
            else
            {
                lnu1 = -0.583 * x;
            }
            nu = Math.Exp(lnu1);
            return nu;
        }

        /// <summary>
        /// integral of 1/(t*(1-t))**2 from x to x+a
        /// </summary>
        /// <param name="x"></param>
        /// <param name="a"></param>
        /// <returns></returns>
        private static double IntegralInvT1tSq(double x, double a)
        {
            double y;
            double integral;

            y = x + a - 0.5;
            integral = (8.0 * y) / (1.0 - 4.0 * Math.Pow(y, 2)) +
                2.0 * Math.Log((1.0 + 2.0 * y) / (1.0 - 2.0 * y));
            y = x - 0.5;
            integral = integral - (8.0 * y) / (1.0 - 4.0 * Math.Pow(y, 2)) -
                2.0 * Math.Log((1.0 + 2.0 * y) / (1.0 - 2.0 * y));
            return integral;
        }


    }
}
