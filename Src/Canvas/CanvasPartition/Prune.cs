using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CanvasPartition
{
    /// <summary>
    /// Contains functions/subroutines found in src/prune.f
    /// </summary>
    internal static class Prune
    {
        /// <summary>
        /// Fortran subroutine errssq
        /// </summary>
        /// <param name="lengthSeg"></param>
        /// <param name="segmentSums">sum of signals for each segment</param>
        /// <param name="k"></param>
        /// <param name="locations"></param>
        /// <returns></returns>
        public static double ErrorSumOfSquares(int[] lengthSeg, double[] segmentSums, int k,
            int[] locations)
        {
            double[] sx = segmentSums; // TODO: remove sx
            double segsx;
            int segnx, i, j;

            double errorSumOfSquares = 0.0;
            segsx = 0.0; // sum of signals across some segments
            segnx = 0; // sum of segment lengths across some segments
            for (i = 0; i < locations[0]; i++)
            {
                segsx += segmentSums[i];
                segnx += lengthSeg[i];
            }
            // (sum of signals across some segments)^2 / (sum of segment lengths across some segments)
            errorSumOfSquares += Math.Pow(segsx, 2) / segnx;
            for (j = 1; j < k; j++)
            {
                segsx = 0.0;
                segnx = 0;
                for (i = locations[j - 1]; i < locations[j]; i++)
                {
                    segsx += segmentSums[i];
                    segnx += lengthSeg[i];
                }
                errorSumOfSquares += Math.Pow(segsx, 2) / segnx;
            }
            segsx = 0.0;
            segnx = 0;
            for (i = locations[k - 1]; i < lengthSeg.Length; i++)
            {
                segsx += segmentSums[i];
                segnx += lengthSeg[i];
            }
            errorSumOfSquares += Math.Pow(segsx, 2) / segnx;
            return errorSumOfSquares;
        }

        /// <summary>
        /// Fortran subroutine combn.
        /// This program generates Choose(n,r) combinations one at a time.
        /// Adapted from Algorithm AS 88  Appl. Statist. (1975) Vol.24, No. 3.
        /// Modifies elements of locations.
        /// Modifies rleft.
        /// </summary>
        /// <param name="r"></param>
        /// <param name="nmr"></param>
        /// <param name="locations"></param>
        /// <param name="rleft"></param>
        public static void Combination(int r, int nmr, int[] locations, ref bool rleft)
        {
            int i, j;

            i = r - 1; // i is an index
            while (locations[i] == nmr + i + 1) { i--; }
            locations[i]++;
            for (j = i + 1; j < r; j++) { locations[j] = locations[j - 1] + 1; }
            if (locations[0] == nmr + 1) { rleft = false; }
        }

    }
}
