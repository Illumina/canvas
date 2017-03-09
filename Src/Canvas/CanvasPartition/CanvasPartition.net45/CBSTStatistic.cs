using System;

namespace CanvasPartition
{
    internal static class CBSTStatistic
    {
        /// <summary>
        /// Fortran subroutine tmaxo
        /// new approach to maximizing t-statistic
        /// dynamic memory allocation using allocatable arrays
        /// </summary>
        /// <param name="genomeData"></param>
        /// <param name="totalSumOfSquares"></param>
        /// <param name="sx"></param>
        /// <param name="iseg"></param>
        /// <param name="ostat"></param>
        /// <param name="al0"></param>
        /// <param name="isBinary"></param>
        public static void TMaxO(double[] genomeData, double totalSumOfSquares, double[] sx,
            int[] iseg, out double ostat, int al0, bool isBinary)
        {
            // look at the partial sums in blocks of size sqrt(n)
            int ipsmin, ipsmax, ipsmin0, ipsmax0, nb, i, j, k, l, nb1,
                nb2, bi, bj, ilo, ihi, jlo, jhi, alenmax, i2j, sxmxi,
                alenlo, alenhi, tmaxi, tmaxj, ixlo, ixhi, nal0;
            double psum, psmin, psmax, psmin0, psmax0, bssmax,
                bsslim, rn, rj, rjhi, rjlo, rnjov1, sij1, sij2, sijmx0,
                absx, sxmx, bijbss, rnov2, psdiff;

            // use local arrays for working within blocks
            // block partial sum max and min 
            double[] bpsmax, bpsmin;
            // location of the max and min
            // bb: block boundary
            int[] bb, ibmin, ibmax;

            // t statistic corresponding to max for block i,j (and max possible)
            // double precision, allocatable :: bssbij(:), bssijmax(:)
            double[] bssbij, bssijmax;
            // row, column and order vector for reordering bssbij
            int[] bloci, blocj, loc, alen;

            // calculate number of blocks (nb) and block boundaries (vector bb)
            rn = Convert.ToDouble(genomeData.Length);
            if (genomeData.Length >= 50)
            {
                nb = Convert.ToInt32(Math.Round(Math.Sqrt(genomeData.Length)));
            }
            else
            {
                nb = 1;
            }
            //the number of pairwise block comparison
            nb2 = nb * (nb + 1) / 2;
            // allocate memory
            bpsmax = new double[nb]; bpsmin = new double[nb];
            bb = new int[nb]; ibmin = new int[nb]; ibmax = new int[nb];
            bssbij = new double[nb2]; bssijmax = new double[nb2];
            bloci = new int[nb2]; blocj = new int[nb2]; loc = new int[nb2]; alen = new int[nb2];

            // block boundaries
            for (i = 0; i < nb; i++)
            {
                bb[i] = Convert.ToInt32(Math.Round(rn * ((i + 1.0) / nb)));
            }
            ilo = 1; // not just an index
            psum = 0;
            psmin0 = 0;
            psmax0 = 0;
            ipsmin0 = genomeData.Length;
            ipsmax0 = genomeData.Length;
            for (j = 0; j < nb; j++) // j is an index and only an index
            {
                sx[ilo - 1] = psum + genomeData[ilo - 1];
                psmin = sx[ilo - 1];
                ipsmin = ilo;
                psmax = sx[ilo - 1];
                ipsmax = ilo;
                for (i = ilo + 1; i <= bb[j]; i++)
                {
                    sx[i - 1] = sx[i - 2] + genomeData[i - 1];
                    if (sx[i - 1] < psmin)
                    {
                        psmin = sx[i - 1];
                        ipsmin = i;
                    }
                    if (sx[i - 1] > psmax)
                    {
                        psmax = sx[i - 1];
                        ipsmax = i;
                    }
                }
                // store the block min, max and locations
                ibmin[j] = ipsmin;
                ibmax[j] = ipsmax;
                bpsmin[j] = psmin;
                bpsmax[j] = psmax;
                // adjust global min, max and locations
                if (psmin < psmin0)
                {
                    psmin0 = psmin;
                    ipsmin0 = ipsmin;
                }
                if (psmax > psmax0)
                {
                    psmax0 = psmax;
                    ipsmax0 = ipsmax;
                }
                // reset ilo to be the block boundary + 1
                psum = sx[bb[j] - 1];
                ilo = bb[j] + 1;
            }
            // calculate bss for max s_i - min s_i
            psdiff = psmax0 - psmin0;
            rj = Math.Abs(ipsmax0 - ipsmin0);
            rnjov1 = rn / (rj * (rn - rj));
            if (isBinary)
            {
                bssmax = rnjov1 * Math.Pow((psdiff - 0.5), 2);
            }
            else
            {
                bssmax = rnjov1 * Math.Pow(psdiff, 2);
            }
            tmaxi = Math.Min(ipsmax0, ipsmin0);
            tmaxj = Math.Max(ipsmax0, ipsmin0);
            // if the segment is all constant then psdiff = 0 and so bssmax = 0
            if (psdiff <= 0)
            {
                bssmax = 0;
                // go to 120
            }
            else
            {
                // for a pair of blocks (i,j) calculate the max absolute t-statistic
                // at the (min_i, max_j) and (max_i, min_j) locations 
                // for other indices the t-statistic can be bounded using this
                //
                //     if a block doesn't have the potential to exceed bssmax ignore it
                //     calculate the bsslim for each block and include ones >= bssmax
                rnov2 = rn / 2;
                l = 0; // l is an index and counter in the following nested for loop
                nal0 = genomeData.Length - al0;
                for (i = 1; i <= nb; i++) // i is not just an 1-based index
                {
                    for (j = i; j <= nb; j++) // j is not just an 1-based index
                    {
                        // calculate bsslim
                        ilo = (i == 1) ? 1 : bb[i - 2] + 1;
                        ihi = bb[i - 1];
                        jlo = (j == 1) ? 1 : bb[j - 2] + 1;
                        jhi = bb[j - 1];
                        alenhi = jhi - ilo;
                        if (alenhi > nal0) { alenhi = nal0; }
                        rjhi = Convert.ToDouble(alenhi);
                        alenlo = (i == j) ? 1 : jlo - ihi;
                        if (alenlo < al0) { alenlo = al0; }
                        // max S_k over block j - min S_k over block i
                        sij1 = Math.Abs(bpsmax[j - 1] - bpsmin[i - 1]);
                        // max S_k over block i - min S_k over block j
                        sij2 = Math.Abs(bpsmax[i - 1] - bpsmin[j - 1]);
                        // if i = j then sij1 and sij2 are the same
                        sijmx0 = Math.Max(sij1, sij2);
                        rjlo = Convert.ToDouble(alenlo);
                        rnjov1 = rn / Math.Min(rjlo * (rn - rjlo), rjhi * (rn - rjhi));
                        if (isBinary)
                        {
                            bsslim = rnjov1 * Math.Pow(sijmx0 - 0.5, 2);
                        }
                        else
                        {
                            bsslim = rnjov1 * Math.Pow(sijmx0, 2);
                        }
                        // if its as large as bssmax add block
                        if (bssmax <= bsslim)
                        {
                            loc[l] = l + 1;
                            bloci[l] = i;
                            blocj[l] = j;
                            bssijmax[l] = bsslim;
                            // max sij in the (i,j) block, t-statistic etc
                            if (sij1 > sij2)
                            {
                                alen[l] = Math.Abs(ibmax[j - 1] - ibmin[i - 1]);
                                rj = Convert.ToDouble(alen[l]);
                                rnjov1 = rn / (rj * (rn - rj));
                                if (isBinary)
                                {
                                    bssbij[l] = rnjov1 * Math.Pow(sij1 - 0.5, 2);
                                }
                                else
                                {
                                    bssbij[l] = rnjov1 * Math.Pow(sij1, 2);
                                }
                            }
                            else
                            {
                                alen[l] = Math.Abs(ibmin[j - 1] - ibmax[i - 1]);
                                rj = Convert.ToDouble(alen[l]);
                                rnjov1 = rn / (rj * (rn - rj));
                                if (isBinary)
                                {
                                    bssbij[l] = rnjov1 * Math.Pow(sij2 - 0.5, 2);
                                }
                                else
                                {
                                    bssbij[l] = rnjov1 * Math.Pow(sij2, 2);
                                }
                            }
                            l++;
                        }
                    }
                }
                nb1 = l;
                // Now sort the t-statistics by their magnitude
                Helper.QuickSort(bssbij, loc, 0, nb1);
                // now go through the blocks in reverse order (largest down)
                for (l = nb1 - 1; l >= 0; l--) // l is an index
                {
                    k = loc[l] - 1; // k is an index in the for loop
                    // need to check a block only if it has potential to increase bss
                    // rjlo is the smalllest (j-i) in the block and rjhi is the largest
                    bsslim = bssijmax[k];
                    if (bssmax <= bsslim)
                    {
                        // bi, bj give the block location
                        bi = bloci[k];
                        bj = blocj[k];
                        // max arc length of interest in block
                        alenmax = alen[k];
                        ilo = (bi == 1) ? 1 : bb[bi - 2] + 1;
                        ihi = bb[bi - 1];
                        jlo = (bj == 1) ? 1 : bb[bj - 2] + 1;
                        jhi = bb[bj - 1];
                        alenhi = jhi - ilo;
                        if (alenhi > nal0) { alenhi = nal0; }
                        rjhi = Convert.ToDouble(alenhi);
                        alenlo = (bi == bj) ? 1 : (jlo - ihi);
                        if (alenlo < al0) { alenlo = al0; }
                        rjlo = Convert.ToDouble(alenlo);
                        // if arc length is larger than n/2 make it n - arc length
                        if (alenmax > genomeData.Length - alenmax)
                        {
                            alenmax = genomeData.Length - alenmax;
                        }
                        // if alenlo <= n/2 start from (ihi, jlo) and go up
                        // if alenhi >= n/2 start from (ilo, jhi) and go down
                        if ((rjlo <= rnov2) && (alenlo <= alenmax))
                        {
                            for (i2j = alenlo; i2j <= alenmax; i2j++)
                            {
                                // excess calculations to set range of i
                                ixlo = Math.Max(0, jlo - ilo - i2j);
                                ixhi = Math.Max(0, ihi + i2j - jhi);
                                sxmx = 0;
                                sxmxi = ilo + ixlo - 1;
                                for (i = ilo + ixlo; i <= ihi - ixhi; i++)
                                {
                                    j = i + i2j;
                                    absx = Math.Abs(sx[j - 1] - sx[i - 1]);
                                    if (sxmx < absx)
                                    {
                                        sxmx = absx;
                                        sxmxi = i;
                                    }
                                }
                                rj = Convert.ToDouble(i2j);
                                rnjov1 = rn / (rj * (rn - rj));
                                if (isBinary)
                                {
                                    bijbss = rnjov1 * Math.Pow(sxmx - 0.5, 2);
                                }
                                else
                                {
                                    bijbss = rnjov1 * Math.Pow(sxmx, 2);
                                }
                                if (bijbss > bssmax)
                                {
                                    bssmax = bijbss;
                                    tmaxi = sxmxi;
                                    tmaxj = sxmxi + i2j;
                                }
                            }
                        }
                        // make arclength n - arc length
                        alenmax = genomeData.Length - alenmax;
                        if ((rjhi >= rnov2) && (alenhi >= alenmax))
                        {
                            for (i2j = alenhi; i2j >= alenmax; i2j--)
                            {
                                // excess calcultaions to set range of i
                                ixlo = Math.Max(0, jlo - ilo - i2j);
                                ixhi = Math.Max(0, ihi + i2j - jhi);
                                sxmx = 0;
                                sxmxi = ilo + ixlo - 1;
                                for (i = ilo + ixlo; i <= ihi - ixhi; i++)
                                {
                                    j = i + i2j;
                                    absx = Math.Abs(sx[j - 1] - sx[i - 1]);
                                    if (sxmx < absx)
                                    {
                                        sxmx = absx;
                                        sxmxi = i;
                                    }
                                }
                                rj = Convert.ToDouble(i2j);
                                rnjov1 = rn / (rj * (rn - rj));
                                if (isBinary)
                                {
                                    bijbss = rnjov1 * Math.Pow(sxmx - 0.5, 2);
                                }
                                else
                                {
                                    bijbss = rnjov1 * Math.Pow(sxmx, 2);
                                }
                                if (bijbss > bssmax)
                                {
                                    bssmax = bijbss;
                                    tmaxi = sxmxi;
                                    tmaxj = sxmxi + i2j;
                                }
                            }
                        }
                    }
                }
            }
            // 120
            if (isBinary)
            {
                if (totalSumOfSquares <= 0.0001) { totalSumOfSquares = 1.0; }
                bssmax = bssmax / (totalSumOfSquares / rn);
            }
            else
            {
                if (totalSumOfSquares <= bssmax + 0.0001) { totalSumOfSquares = bssmax + 1.0; }
                bssmax = bssmax / ((totalSumOfSquares - bssmax) / (rn - 2.0));
            }
            ostat = bssmax;
            iseg[0] = tmaxi;
            iseg[1] = tmaxj;
        }

        /// <summary>
        /// function for the max (over small arcs) t-statistic on permuted data
        /// new code to speed up this part 3/31/2010
        /// </summary>
        /// <param name="k"></param>
        /// <param name="totalSumOfSquares">total sum of squares</param>
        /// <param name="px">permuted genome data</param>
        /// <param name="sx"></param>
        /// <param name="al0"></param>
        /// <param name="isBinary">is the data binary?</param>
        /// <returns></returns>
        public static double HTMaxP(int k, double totalSumOfSquares, double[] px, double[] sx,
            int al0, bool isBinary)
        {
            int i, j, nmj;
            double rn, rj, absx, sxmx, bssmx, psmin, psmax, psdiff,
                bsslim, rnjov1;
            int ipsmin, ipsmax;

            // create blocks of size k (or k+1) to span 1 thru n
            // block partial sum max and min 
            double[] bpsmax, bpsmin;
            // location of the max and min
            int[] bb;
            // variables to work on block specific data
            int nb, ilo, ihi, l;
            double psum, psdiffsq;

            rn = Convert.ToDouble(px.Length);
            // number of blocks of size k (plus fraction since n/k may not be integer)
            nb = (int)(rn / k);
            // allocate memory
            bpsmax = new double[nb]; bpsmin = new double[nb];
            bb = new int[nb];
            // block boundaries
            for (i = 0; i < nb; i++)
            {
                bb[i] = Convert.ToInt32(Math.Round(rn * (Convert.ToDouble(i + 1) / nb)));
            }

            // don't need global min and max
            // find the max, min of partial sums and their locations within blocks
            ilo = 1;
            psum = 0;
            double hTMaxP = 0.0;
            for (j = 0; j < nb; j++) // j is just an index in this for loop
            {
                sx[ilo - 1] = psum + px[ilo - 1];
                psmin = sx[ilo - 1];
                ipsmin = ilo;
                psmax = sx[ilo - 1];
                ipsmax = ilo;
                for (i = ilo; i < bb[j]; i++)
                {
                    sx[i] = sx[i - 1] + px[i];
                    if (sx[i] < psmin)
                    {
                        psmin = sx[i];
                        ipsmin = i + 1;
                    }
                    if (sx[i] > psmax)
                    {
                        psmax = sx[i];
                        ipsmax = i + 1;
                    }
                }
                // store the block min, max and locations
                bpsmin[j] = psmin;
                bpsmax[j] = psmax;
                // reset ilo to be the block boundary + 1
                psum = sx[bb[j] - 1];
                ilo = bb[j] + 1;
                // calculate the bss at the block max & min pr
                i = Math.Abs(ipsmin - ipsmax);
                if ((i <= k) && (i >= al0))
                {
                    rj = Convert.ToDouble(i);
                    rnjov1 = rn / (rj * (rn - rj));
                    if (isBinary)
                    {
                        bssmx = rnjov1 * Math.Pow(bpsmax[j] - bpsmin[j] - 0.5, 2);
                    }
                    else
                    {
                        bssmx = rnjov1 * Math.Pow(bpsmax[j] - bpsmin[j], 2);
                    }
                    if (hTMaxP < bssmx) { hTMaxP = bssmx; }
                }
            }
            // check the first block
            ilo = 1;
            ihi = bb[0];
            psdiff = bpsmax[0] - bpsmin[0];
            if (isBinary)
            {
                psdiffsq = Math.Pow(psdiff - 0.5, 2);
            }
            else
            {
                psdiffsq = Math.Pow(psdiff, 2);
            }
            for (j = al0; j <= k; j++)
            {
                rj = Convert.ToDouble(j);
                rnjov1 = rn / (rj * (rn - rj));
                bsslim = rnjov1 * psdiffsq;
                if (bsslim < hTMaxP) { break; } // go to 50
                sxmx = 0.0;
                for (i = ilo; i <= ihi - j; i++)
                {
                    absx = Math.Abs(sx[i + j - 1] - sx[i - 1]);
                    if (sxmx < absx) { sxmx = absx; }
                }
                if (isBinary)
                {
                    bssmx = rnjov1 * Math.Pow(Math.Abs(sxmx) - 0.5, 2);
                }
                else
                {
                    bssmx = rnjov1 * Math.Pow(sxmx, 2);
                }
                if (hTMaxP < bssmx) { hTMaxP = bssmx; }
            }
            // 50 now the minor arcs spanning the end (n)
            psdiff = Math.Max(Math.Abs(bpsmax[0] - bpsmin[nb - 1]),
                Math.Abs(bpsmax[nb - 1] - bpsmin[0]));
            if (isBinary)
            {
                psdiffsq = Math.Pow(psdiff - 0.5, 2);
            }
            else
            {
                psdiffsq = Math.Pow(psdiff, 2);
            }
            for (j = al0; j <= k; j++)
            {
                rj = Convert.ToDouble(j);
                rnjov1 = rn / (rj * (rn - rj));
                bsslim = rnjov1 * psdiffsq;
                if (bsslim < hTMaxP) { break; } // go to 100
                sxmx = 0.0;
                nmj = px.Length - j;
                for (i = 0; i < j; i++) // i is just an index in this for loop
                {
                    absx = Math.Abs(sx[i + nmj] - sx[i]);
                    if (sxmx < absx) { sxmx = absx; }
                }
                if (isBinary)
                {
                    bssmx = rnjov1 * Math.Pow(Math.Abs(sxmx) - 0.5, 2);
                }
                else
                {
                    bssmx = rnjov1 * Math.Pow(sxmx, 2);
                }
                if (hTMaxP < bssmx) { hTMaxP = bssmx; }
            }
            // 100 now the other blocks
            for (l = 1; l < nb; l++) // l is just an index in this for loop
            {
                ilo = bb[l - 1] + 1;
                ihi = bb[l];
                psdiff = bpsmax[l] - bpsmin[l];
                if (isBinary)
                {
                    psdiffsq = Math.Pow(psdiff - 0.5, 2);
                }
                else
                {
                    psdiffsq = Math.Pow(psdiff, 2);
                }
                for (j = al0; j <= k; j++)
                {
                    rj = Convert.ToDouble(j);
                    rnjov1 = rn / (rj * (rn - rj));
                    bsslim = rnjov1 * psdiffsq;
                    if (bsslim < hTMaxP) { break; } // go to 150
                    sxmx = 0.0;
                    for (i = ilo; i <= ihi - j; i++)
                    {
                        absx = Math.Abs(sx[i + j - 1] - sx[i - 1]);
                        if (sxmx < absx) { sxmx = absx; }
                    }
                    if (isBinary)
                    {
                        bssmx = rnjov1 * Math.Pow(Math.Abs(sxmx) - 0.5, 2);
                    }
                    else
                    {
                        bssmx = rnjov1 * Math.Pow(sxmx, 2);
                    }
                    if (hTMaxP < bssmx) { hTMaxP = bssmx; }
                }
                // 150
                psdiff = Math.Max(Math.Abs(bpsmax[l] - bpsmin[l - 1]),
                    Math.Abs(bpsmax[l - 1] - bpsmin[l]));
                if (isBinary)
                {
                    psdiffsq = Math.Pow(psdiff - 0.5, 2);
                }
                else
                {
                    psdiffsq = Math.Pow(psdiff, 2);
                }
                for (j = al0; j <= k; j++)
                {
                    rj = Convert.ToDouble(j);
                    rnjov1 = rn / (rj * (rn - rj));
                    bsslim = rnjov1 * psdiffsq;
                    if (bsslim < hTMaxP) { break; } // go to 200
                    sxmx = 0.0;
                    nmj = px.Length - j;
                    for (i = ilo - j; i <= ilo - 1; i++)
                    {
                        absx = Math.Abs(sx[i + j - 1] - sx[i - 1]);
                        if (sxmx < absx) { sxmx = absx; }
                    }
                    if (isBinary)
                    {
                        bssmx = rnjov1 * Math.Pow(Math.Abs(sxmx) - 0.5, 2);
                    }
                    else
                    {
                        bssmx = rnjov1 * Math.Pow(sxmx, 2);
                    }
                    if (hTMaxP < bssmx) { hTMaxP = bssmx; }
                }
            }// 200
            if (isBinary)
            {
                if (totalSumOfSquares <= 0.0001) { totalSumOfSquares = 1.0; }
                hTMaxP = hTMaxP / (totalSumOfSquares / rn);
            }
            else
            {
                if (totalSumOfSquares <= hTMaxP + 0.0001)
                {
                    totalSumOfSquares = hTMaxP + 1.0;
                }
                hTMaxP = hTMaxP / ((totalSumOfSquares - hTMaxP) / (rn - 2.0));
            }

            return hTMaxP;
        }


        /// <summary>
        /// function for calculating the full max t-statistic on permuted data
        /// new approach to maximizing t-statistic using allocatable arrays 
        /// </summary>
        /// <param name="totalSumOfSquares"></param>
        /// <param name="px">permuted genome data</param>
        /// <param name="sx"></param>
        /// <param name="al0"></param>
        /// <param name="isBinary"></param>
        /// <returns></returns>
        public static double TMaxP(double totalSumOfSquares, double[] px, double[] sx,
            int al0, bool isBinary)
        {
            // look at the partial sums in blocks of size sqrt(n)

            int ipsmin, ipsmax, ipsmin0, ipsmax0, nb, i, j, k, l, nb1,
                nb2, bi, bj, ilo, ihi, jlo, jhi, alenmax, i2j, alenlo,
                alenhi, ixlo, ixhi, nal0;
            double psum, psmin, psmax, psmin0, psmax0, bssmax,
                bsslim, rn, rj, rjhi, rjlo, rnjov1, sij1, sij2, sijmx0,
                absx, sxmx, bijbss, rnov2, psdiff;

            // use local arrays for working within blocks
            // block partial sum max and min 
            double[] bpsmax, bpsmin;
            // location of the max and min
            int[] bb, ibmin, ibmax;

            // t statistic corresponding to max for block i,j (and max possible)
            double[] bssbij, bssijmax;
            // row, column and order vector for reordering bssbij
            int[] bloci, blocj, loc, alen;

            // calculate number of blocks (nb) and block boundaries (vector bb)
            rn = Convert.ToDouble(px.Length);
            if (px.Length >= 50)
            {
                nb = Convert.ToInt32(Math.Round(Math.Sqrt(px.Length)));
            }
            else
            {
                nb = 1;
            }

            // the number of paiwise block comparison
            nb2 = nb * (nb + 1) / 2;
            // allocate memory
            bpsmax = new double[nb]; bpsmin = new double[nb];
            bb = new int[nb]; ibmin = new int[nb]; ibmax = new int[nb];
            bssbij = new double[nb2]; bssijmax = new double[nb2];
            bloci = new int[nb2]; blocj = new int[nb2]; loc = new int[nb2]; alen = new int[nb2];

            // block boundaries
            for (i = 0; i < nb; i++)
            {
                bb[i] = Convert.ToInt32(Math.Round(rn * (Convert.ToDouble(i + 1) / nb)));
            }

            // find the max, min of partial sums and their locations within blocks
            ilo = 1;
            psum = 0;
            psmin0 = 0;
            psmax0 = 0;
            ipsmin0 = px.Length;
            ipsmax0 = px.Length;
            for (j = 0; j < nb; j++) // j is just an index
            {
                sx[ilo - 1] = psum + px[ilo - 1];
                psmin = sx[ilo - 1];
                ipsmin = ilo;
                psmax = sx[ilo - 1];
                ipsmax = ilo;
                for (i = ilo; i < bb[j]; i++)
                {
                    sx[i] = sx[i - 1] + px[i];
                    if (sx[i] < psmin)
                    {
                        psmin = sx[i];
                        ipsmin = i + 1;
                    }
                    if (sx[i] > psmax)
                    {
                        psmax = sx[i];
                        ipsmax = i + 1;
                    }
                }
                // store the block min, max and locations
                ibmin[j] = ipsmin;
                ibmax[j] = ipsmax;
                bpsmin[j] = psmin;
                bpsmax[j] = psmax;
                // adjust global min, max and locations
                if (psmin < psmin0)
                {
                    psmin0 = psmin;
                    ipsmin0 = ipsmin;
                }
                if (psmax > psmax0)
                {
                    psmax0 = psmax;
                    ipsmax0 = ipsmax;
                }
                // reset ilo to be the block boundary + 1
                psum = sx[bb[j] - 1];
                ilo = bb[j] + 1;
            }

            // calculate bss for max s_i - min s_i
            psdiff = psmax0 - psmin0;
            rj = Math.Abs(ipsmax0 - ipsmin0);
            rnjov1 = rn / (rj * (rn - rj));
            if (isBinary)
            {
                bssmax = rnjov1 * Math.Pow(psdiff - 0.5, 2);
            }
            else
            {
                bssmax = rnjov1 * Math.Pow(psdiff, 2);
            }

            // for a pair of blocks (i,j) calculate the max absolute t-statistic
            // at the (min_i, max_j) and (max_i, min_j) locations 
            // for other indices the t-statistic can be bounded using this

            // if a block doesn't have the potential to exceed bssmax ignore it
            // calculate the bsslim for each block and include ones >= bssmax

            rnov2 = rn / 2;
            l = 0;
            nal0 = px.Length - al0;
            for (i = 1; i <= nb; i++)
            {
                for (j = i; j <= nb; j++)
                {
                    // calculate bsslim
                    if (i == 1)
                    {
                        ilo = 1;
                    }
                    else
                    {
                        ilo = bb[i - 2] + 1;
                    }
                    ihi = bb[i - 1];
                    if (j == 1)
                    {
                        jlo = 1;
                    }
                    else
                    {
                        jlo = bb[j - 2] + 1;
                    }
                    jhi = bb[j - 1];
                    alenhi = jhi - ilo;
                    if (alenhi > nal0) { alenhi = nal0; }
                    rjhi = Convert.ToDouble(alenhi);
                    if (i == j)
                    {
                        alenlo = 1;
                    }
                    else
                    {
                        alenlo = jlo - ihi;
                    }
                    if (alenlo < al0) { alenlo = al0; }
                    // max S_k over block j - min S_k over block i
                    sij1 = Math.Abs(bpsmax[j - 1] - bpsmin[i - 1]);
                    // max S_k over block i - min S_k over block j
                    sij2 = Math.Abs(bpsmax[i - 1] - bpsmin[j - 1]);
                    // if i = j then sij1 and sij2 are the same
                    sijmx0 = Math.Max(sij1, sij2);
                    rjlo = Convert.ToDouble(alenlo);
                    rnjov1 = rn / Math.Min(rjlo * (rn - rjlo), rjhi * (rn - rjhi));
                    if (isBinary)
                    {
                        bsslim = rnjov1 * Math.Pow(sijmx0 - 0.5, 2);
                    }
                    else
                    {
                        bsslim = rnjov1 * Math.Pow(sijmx0, 2);
                    }
                    // if its as large as bssmax add block
                    if (bssmax <= bsslim)
                    {
                        loc[l] = l;
                        bloci[l] = i;
                        blocj[l] = j;
                        bssijmax[l] = bsslim;
                        // max sij in the (i,j) block, t-statistic etc
                        if (sij1 > sij2)
                        {
                            alen[l] = Math.Abs(ibmax[j - 1] - ibmin[i - 1]);
                            rj = Convert.ToDouble(alen[l]);
                            rnjov1 = rn / (rj * (rn - rj));
                            if (isBinary)
                            {
                                bssbij[l] = rnjov1 * Math.Pow(sij1 - 0.5, 2);
                            }
                            else
                            {
                                bssbij[l] = rnjov1 * Math.Pow(sij1, 2);
                            }
                        }
                        else
                        {
                            alen[l] = Math.Abs(ibmin[j - 1] - ibmax[i - 1]);
                            rj = Convert.ToDouble(alen[l]);
                            rnjov1 = rn / (rj * (rn - rj));
                            if (isBinary)
                            {
                                bssbij[l] = rnjov1 * Math.Pow(sij2 - 0.5, 2);
                            }
                            else
                            {
                                bssbij[l] = rnjov1 * Math.Pow(sij2, 2);
                            }
                        }
                        l = l + 1;
                    }
                }
            }
            nb1 = l;
            // Now sort the t-statistics by their magnitude
            Helper.QuickSort(bssbij, loc, 0, nb1);

            // now go through the blocks in reverse order (largest down)
            for (l = nb1 - 1; l >= 0; l--) // l is an index in this for loop
            {
                k = loc[l] - 1; // k is an index in this for loop
                // need to check a block only if it has potential to increase bss
                // rjlo is the smalllest (j-i) in the block and rjhi is the largest
                bsslim = bssijmax[k];
                if (bssmax <= bsslim)
                {
                    // bi, bj give the block location
                    bi = bloci[k];
                    bj = blocj[k];
                    // max arc length of interest in block
                    alenmax = alen[k];
                    if (bi == 1)
                    {
                        ilo = 1;
                    }
                    else
                    {
                        ilo = bb[bi - 2] + 1;
                    }
                    ihi = bb[bi - 1];
                    if (bj == 1)
                    {
                        jlo = 1;
                    }
                    else
                    {
                        jlo = bb[bj - 2] + 1;
                    }
                    jhi = bb[bj - 1];
                    alenhi = jhi - ilo;
                    if (alenhi > nal0) { alenhi = nal0; }
                    rjhi = Convert.ToDouble(alenhi);
                    if (bi == bj)
                    {
                        alenlo = 1;
                    }
                    else
                    {
                        alenlo = jlo - ihi;
                    }
                    if (alenlo < al0) { alenlo = al0; }
                    rjlo = Convert.ToDouble(alenlo);
                    // if arc length is larger than n/2 make is n - arc length
                    if (alenmax > px.Length - alenmax) { alenmax = px.Length - alenmax; }
                    // if alenlo <= n/2 start from (ihi, jlo) and go up
                    // if alenhi >= n/2 start from (ilo, jhi) and go down
                    if ((rjlo <= rnov2) && (alenlo <= alenmax))
                    {
                        for (i2j = alenlo; i2j <= alenmax; i2j++)
                        {
                            // excess calcultaions to set range of i
                            ixlo = Math.Max(0, jlo - ilo - i2j);
                            ixhi = Math.Max(0, ihi + i2j - jhi);
                            sxmx = 0;
                            for (i = ilo + ixlo - 1; i < ihi - ixhi; i++) //i: 0-based index
                            {
                                j = i + i2j; // j: 0-based index
                                absx = Math.Abs(sx[j] - sx[i]);
                                if (sxmx < absx) { sxmx = absx; }
                            }
                            rj = Convert.ToDouble(i2j);
                            rnjov1 = rn / (rj * (rn - rj));
                            if (isBinary)
                            {
                                bijbss = rnjov1 * Math.Pow(sxmx - 0.5, 2);
                            }
                            else
                            {
                                bijbss = rnjov1 * Math.Pow(sxmx, 2);
                            }
                            if (bijbss > bssmax) { bssmax = bijbss; }
                        }
                    }
                    // make arclength n - arc length
                    alenmax = px.Length - alenmax;
                    if ((rjhi >= rnov2) && (alenhi >= alenmax))
                    {
                        for (i2j = alenhi; i2j >= alenmax; i2j--)
                        {
                            // excess calcultaions to set range of i
                            ixlo = Math.Max(0, jlo - ilo - i2j);
                            ixhi = Math.Max(0, ihi + i2j - jhi);
                            sxmx = 0;
                            for (i = ilo + ixlo - 1; i < ihi - ixhi; i++)  // i: 0-based index
                            {
                                j = i + i2j; // j: 0-based index
                                absx = Math.Abs(sx[j] - sx[i]);
                                if (sxmx < absx) { sxmx = absx; }
                            }
                            rj = Convert.ToDouble(i2j);
                            rnjov1 = rn / (rj * (rn - rj));
                            if (isBinary)
                            {
                                bijbss = rnjov1 * Math.Pow(sxmx - 0.5, 2);
                            }
                            else
                            {
                                bijbss = rnjov1 * Math.Pow(sxmx, 2);
                            }
                            if (bijbss > bssmax) { bssmax = bijbss; }
                        }
                    }
                }
            }
            double tMaxP = 0.0;
            if (isBinary)
            {
                if (totalSumOfSquares <= 0.0001) { totalSumOfSquares = 1.0; }
                tMaxP = bssmax / (totalSumOfSquares / rn);
            }
            else
            {
                if (totalSumOfSquares <= bssmax + 0.0001) { totalSumOfSquares = bssmax + 1.0; }
                tMaxP = bssmax / ((totalSumOfSquares - bssmax) / (rn - 2.0));
            }

            return tMaxP;
        }


        /// <summary>
        /// function for the p-value of t-statistics for removing edge effects
        /// </summary>
        /// <param name="n1"></param>
        /// <param name="n2"></param>
        /// <param name="genomeData"></param>
        /// <param name="gDOffset">starts from genomeData[gDOffset]</param>
        /// <param name="px">permuted genome data</param>
        /// <param name="nPerm">number of permutation</param>
        /// <returns></returns>
        public static double TPermP(int n1, int n2, int n, double[] genomeData, int gDOffset,
            double[] px, uint nPerm, Random rnd)
        {
            int np, i, m1, j, nrej;
            double xsum1, xsum2, xbar, ostat, pstat, rn1, rn2, rm1,
                tstat, tss, rn, cc;

            rn1 = Convert.ToDouble(n1);
            rn2 = Convert.ToDouble(n2);
            rn = rn1 + rn2;
            if (n1 == 1 || n2 == 1)
            {
                nrej = Convert.ToInt32(nPerm);
                // go to 110
            }
            else
            {
                xsum1 = 0.0;
                tss = 0.0;
                for (i = 0; i < n1; i++) // i: 0-based index in the for loop
                {
                    px[i] = genomeData[gDOffset + i];
                    xsum1 = xsum1 + genomeData[gDOffset + i];
                    tss = tss + Math.Pow(genomeData[gDOffset + i], 2);
                }
                xsum2 = 0.0;
                for (i = n1; i < n; i++)
                {
                    px[i] = genomeData[gDOffset + i];
                    xsum2 = xsum2 + genomeData[gDOffset + i];
                    tss = tss + Math.Pow(genomeData[gDOffset + i], 2);
                }
                xbar = (xsum1 + xsum2) / rn;
                tss = tss - rn * Math.Pow(xbar, 2);
                if (n1 <= n2)
                {
                    m1 = n1;
                    rm1 = rn1;
                    ostat = 0.99999 * Math.Abs(xsum1 / rn1 - xbar);
                    tstat = Math.Pow(ostat, 2) * rn1 * rn / rn2;
                }
                else
                {
                    m1 = n2;
                    rm1 = rn2;
                    ostat = 0.99999 * Math.Abs(xsum2 / rn2 - xbar);
                    tstat = Math.Pow(ostat, 2) * rn2 * rn / rn1;
                }
                nrej = 0;
                tstat = tstat / ((tss - tstat) / (rn - 2.0));
                // if observed t is large (> 5) don't bother with permutation p-value
                // also make sure there are enough observations i.e. m1 >= 10
                if ((tstat > 25) && (m1 >= 10))
                { // go to 110
                }
                else
                {
                    for (np = 0; np < nPerm; np++)
                    {
                        xsum1 = 0;
                        for (i = n - 1; i >= n - m1; i--)
                        {
                            cc = rnd.NextDouble();
                            j = (int)(cc * (i + 1));
                            j = (j > i) ? i : j;
                            Helper.Swap<double>(ref px[i], ref px[j]);
                            xsum1 = xsum1 + px[i];
                        }
                        pstat = Math.Abs(xsum1 / rm1 - xbar);
                        if (ostat <= pstat) { nrej = nrej + 1; }
                    }
                }
            }
            // 110 
            double tPermP = Convert.ToDouble(nrej) / nPerm;

            return tPermP;
        }




    }
}
