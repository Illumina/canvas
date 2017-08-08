using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.Distributions;

namespace CanvasPartition
{
    internal enum SegmentSplitUndo
    {
        None = 0,
        Prune,
        SDUndo
    }


    /// <summary>
    /// Contains static funtions found in R/changepoints.R and src/changepoints.f
    /// </summary>
    internal static class ChangePoint
    {
        /// <summary>
        /// Outputs:
        ///     lengthSeg
        ///     segmentMeans
        /// </summary>
        /// <param name="genomeData"></param>
        /// <param name="sbdry"></param>
        /// <param name="lengthSeg">segment lengths</param>
        /// <param name="segmentMeans">segment means</param>
        /// <param name="dataType">"logratio" or "binary"</param>
        /// <param name="alpha"></param>
        /// <param name="nPerm"></param>
        /// <param name="pMethod"></param>
        /// <param name="minWidth"></param>
        /// <param name="kMax"></param>
        /// <param name="nMin"></param>
        /// <param name="trimmedSD"></param>
        /// <param name="undoSplits">"none" or "prune" or "sdundo"</param>
        /// <param name="undoPrune"></param>
        /// <param name="undoSD"></param>
        /// <param name="verbose"></param>
        /// <param name="nGrid"></param>
        /// <param name="tol"></param>
        public static void ChangePoints(double[] genomeData, uint[] sbdry,
            out int[] lengthSeg, out double[] segmentMeans, Random rnd,
            string dataType = "logratio", double alpha = 0.01, 
            uint nPerm = 10000, string pMethod = "hybrid", int minWidth = 2, int kMax = 25,
            uint nMin = 200, double trimmedSD = -1, SegmentSplitUndo undoSplits = SegmentSplitUndo.None,
            double undoPrune = 0.05, double undoSD = 3, int verbose = 1,
            int nGrid = 100, double tol = 1E-6)
        {
            if (trimmedSD <= 0)
            {
                trimmedSD = Helper.MedianAbsoluteDeviation(Helper.Diff(genomeData)) / Math.Sqrt(2);
            }
            // start with the whole 
            List<int> segEnd = new List<int>();
            segEnd.Add(0); // inclusive
            segEnd.Add(genomeData.Length); // exclusive
            int k = segEnd.Count;
            List<int> changeLocations = new List<int>();
            int nChangePoints = 0;
            int[] iChangePoint = null;
            while (k > 1)
            {
                int currentN = segEnd[k - 1] - segEnd[k - 2];
                if (verbose >= 3)
                {
                    Console.Write(".... current segment: {0} - {1} \n", segEnd[k - 2] + 1, segEnd[k - 1]);
                }
                if (currentN >= 2 * minWidth)
                {
                    double[] currentGenomeData = new double[currentN];
                    Array.Copy(genomeData, segEnd[k - 2], currentGenomeData, 0, currentN);
                    // check whether hybrid method needs to be used
                    bool hybrid = false;
                    double delta = 0.0;
                    if (pMethod.Equals("hybrid") && nMin < currentN)
                    {
                        hybrid = true;
                        delta = (kMax + 1.0) / currentN;
                    }

                    // if all values of current.genomdat are the same don't segment
                    if (currentGenomeData.Max() == currentGenomeData.Min())
                    {
                        nChangePoints = 0;
                    }
                    else
                    {
                        // centering the current data will save a lot of computations later
                        double currentAverage = currentGenomeData.Average();
                        Helper.InplaceSub(currentGenomeData, currentAverage);
                        // need total sum of squares too
                        double currentTSS = Helper.WeightedSumOfSquares(currentGenomeData, null);
                        ChangePoint.FindChangePoints(currentGenomeData, currentTSS, nPerm, alpha,
                            out nChangePoints, out iChangePoint, dataType.Equals("binary"),
                            hybrid, minWidth, kMax, delta, nGrid, sbdry, tol, rnd);
                    }
                }
                else
                {
                    nChangePoints = 0;
                }
                // Save the change location
                // segEnd[k - 1] will be removed when nChangePoints == 0
                if (nChangePoints == 0) { changeLocations.Add(segEnd[k - 1]); }
                // Offset iChangePoint by segEnd[k - 2]
                for (int i = 0; i < nChangePoints; i++) { iChangePoint[i] += segEnd[k - 2]; }
                switch (nChangePoints) // switch by the number of change points
                {
                    case 0: // no change point
                        segEnd.RemoveAt(k - 1); // Remove the last element
                        break;
                    case 1: // one change point
                        segEnd.Insert(k - 1, iChangePoint[0]);
                        break;
                    case 2: // two change points
                        segEnd.InsertRange(k - 1, iChangePoint);
                        break;
                    default:
                        Console.Error.WriteLine("There should be 0, 1, or 2 change points");
                        break;
                }
                k = segEnd.Count;
                if (verbose >= 3) { Console.Write(".... segments to go: {0} \n", String.Join(" ", segEnd)); }
            }
            changeLocations.Reverse(); // changeLocations is no longer needed
            List<int> segEnds = changeLocations;
            int nSeg = segEnds.Count;
            segEnds.Insert(0, 0);
            lengthSeg = Helper.Diff(segEnds.ToArray());
            if (nSeg > 1)
            {
                if (undoSplits == SegmentSplitUndo.Prune)
                {
                    lengthSeg = ChangePointsPrune(genomeData, lengthSeg, changeCutoff: undoPrune);
                }
                if (undoSplits == SegmentSplitUndo.SDUndo)
                {
                    lengthSeg = ChangePointsSDUndo(genomeData, lengthSeg, trimmedSD, changeSD: undoSD);
                }
            }
            segmentMeans = new double[lengthSeg.Length];
            int ll = 0, uu = 0;
            for (int i = 0; i < lengthSeg.Length; i++)
            {
                uu += lengthSeg[i];
                // Works even if weights == null
                segmentMeans[i] = Helper.WeightedAverage(genomeData, null, iStart: ll, iEnd: uu);
                ll = uu;
            }
        }

        private static int[] ChangePointsSDUndo(double[] genomeData, int[] lengthSeg, double trimmedSD,
            double changeSD = 3)
        {
            changeSD *= trimmedSD;
            var changePointLocations = new List<int>(Helper.CumulativeSum(lengthSeg));
            bool sdUndo = true;
            while (sdUndo)
            {
                int k = changePointLocations.Count;
                if (k > 1)
                {
                    var starts = new List<int>(changePointLocations); // make a copy of changePointLocations
                    starts.RemoveAt(k - 1); // Remove the element at k - 1
                    starts.Insert(0, 0); // Insert 0 to the front of tmp
                    // starts will be used as the start indices (0-based, inclusive) ==> no need to add 1
                    // changePointLocations will be used as the end indices (0-based, exclusive)
                    double[] segmentMedians = new double[k];
                    for (int i = 0; i < k; i++) // for each segment
                    {
                        segmentMedians[i] = Helper.Median(genomeData, iStart: starts[i],
                            iEnd: changePointLocations[i]);
                    }
                    double[] absDiffSegmentMedians = Helper.InplaceAbs(Helper.Diff(segmentMedians));
                    double min;
                    int iMin = Helper.ArgMin(absDiffSegmentMedians, out min);
                    if (min < changeSD)
                    {
                        changePointLocations.RemoveAt(iMin);
                    }
                    else
                    {
                        sdUndo = false;
                    }
                }
                else
                {
                    sdUndo = false;
                }
            }
            changePointLocations.Insert(0, 0);
            return Helper.Diff(changePointLocations.ToArray()); // segment lengths
        }

        /// <summary>
        /// R function changepoints.prune, which calls Fortran subroutine prune
        /// </summary>
        /// <param name="genomeData"></param>
        /// <param name="lengthSeg"></param>
        /// <param name="changeCutoff"></param>
        /// <returns></returns>
        private static int[] ChangePointsPrune(double[] genomeData, int[] lengthSeg, double changeCutoff = 0.05)
        {
            double[] sx = new double[lengthSeg.Length]; // segment sums
            int[] loc = new int[lengthSeg.Length - 1]; // one element for each change point
            int[,] loc1 = new int[2, lengthSeg.Length - 1]; // one element for each change point
            int prunedNChangePoints = 0; // number of change points after pruning

            int i, j, k, kmj; // kmj: k minus j
            double ssq, wssqk, wssq1, wssqj;
            bool jleft;

            ssq = Helper.PartialSumOfPowers(genomeData, 2, 0, genomeData.Length); // sum of squares
            k = 0; // segment start index 
            for (i = 0; i < lengthSeg.Length; i++) // for each segment
            {
                sx[i] = Helper.PartialSumOfPowers(genomeData, 1, k, lengthSeg[i]);
                k += lengthSeg[i];
            }
            // k = nseg - 1 == number of change points
            for (i = 0; i < loc.Length; i++)
            {
                loc[i] = i + 1;
                loc1[1, i] = i + 1;
            }
            wssqk = ssq - Prune.ErrorSumOfSquares(lengthSeg, sx, loc.Length, loc);
            // j: number of change points
            for (j = loc.Length - 1; j > 0; j--) // j (= loc.Length - 1, ..., 1) is not an index
            {
                kmj = loc.Length - j;
                jleft = true;
                for (i = 0; i < j; i++)
                {
                    loc[i] = i + 1;
                    loc1[0, i] = i + 1;
                }
                wssqj = ssq - Prune.ErrorSumOfSquares(lengthSeg, sx, j, loc);
                while (jleft)
                {
                    Prune.Combination(j, kmj, loc, ref jleft);
                    wssq1 = ssq - Prune.ErrorSumOfSquares(lengthSeg, sx, j, loc);
                    if (wssq1 <= wssqj)
                    {
                        wssqj = wssq1;
                        for (i = 0; i < j; i++) { loc1[0, i] = loc[i]; }
                    }
                }
                if (wssqj / wssqk > 1 + changeCutoff)
                {
                    prunedNChangePoints = j + 1; // go back to the previous j
                    for (i = 0; i < prunedNChangePoints; i++) { loc[i] = loc1[1, i]; }
                    break;
                }
                else
                {
                    for (i = 0; i < j; i++) { loc1[1, i] = loc1[0, i]; }
                }
            }
            int[] cumSumLengthSeg = Helper.CumulativeSum(lengthSeg);
            int[] prunedChangePoints = new int[prunedNChangePoints + 2];
            for (i = 0; i < prunedNChangePoints; i++)
            {   // loc[i] is an 1-based index
                prunedChangePoints[i + 1] = cumSumLengthSeg[loc[i] - 1];
            }
            prunedChangePoints[0] = 0;
            prunedChangePoints[prunedChangePoints.Length - 1] = genomeData.Length;
            return Helper.Diff(prunedChangePoints);
        }

        /// <summary>
        /// Fortran subroutine fndcpt.
        /// Ternary segmentation with permutation reference distribution
        /// </summary>
        /// <param name="genomeData"></param>
        /// <param name="totalSumOfSquares"></param>
        /// <param name="nPerm"></param>
        /// <param name="cutoffPValue"></param>
        /// <param name="nChangePoints"></param>
        /// <param name="iChangePoint"></param>
        /// <param name="isBinary"></param>
        /// <param name="hybrid"></param>
        /// <param name="al0"></param>
        /// <param name="hk"></param>
        /// <param name="delta"></param>
        /// <param name="nGrid"></param>
        /// <param name="sbdry"></param>
        /// <param name="tol"></param>
        private static void FindChangePoints(double[] genomeData, double totalSumOfSquares, uint nPerm,
            double cutoffPValue, out int nChangePoints, out int[] iChangePoint, bool isBinary,
            bool hybrid, int al0, int hk, double delta, int nGrid, uint[] sbdry, double tol, Random rnd)
        {
            double[] px = new double[genomeData.Length]; // permuted genomeData
            double[] sx = new double[genomeData.Length];
            iChangePoint = new int[2]; // up to 2 change points

            // nrej: # of non-rejected tests
            // nrejc: # of non-rejected tests cutoff
            int np, nrej, nrejc, n1, n2, n12, l, k;
            int[] iseg = new int[2];
            // segment lengths: iseg[0], iseg[1] - iseg[0], genomeData.Length - iseg[1]
            double ostat, ostat1, pstat, tPValue, pValue1, pValue2;

            nrej = 0;
            nChangePoints = 0;

            CBSTStatistic.TMaxO(genomeData, totalSumOfSquares, sx, iseg, out ostat, al0, isBinary);
            ostat1 = Math.Sqrt(ostat);
            ostat *= 0.99999;
            // if maximal t-statistic is too small (for now use 0.1) don't split
            if (ostat1 <= 0.1) { return; } // call rndend() before return?
            // if maximal t-statistic is too large (for now use 7.0) split
            // also make sure it's not affected by outliers i.e. small seglength
            l = Math.Min(iseg[1] - iseg[0], genomeData.Length - iseg[1] + iseg[0]);
            if (!((ostat1 >= 7.0) && (l >= 10)))
            {
                // o.w calculate p-value and decide if & how data are segmented
                if (hybrid)
                {
                    pValue1 = TailProbability.TailP(ostat1, delta, genomeData.Length, nGrid, tol);
                    if (pValue1 > cutoffPValue) { return; } // pValue1 is the lower bound
                    pValue2 = cutoffPValue - pValue1;
                    nrejc = (int)(pValue2 * nPerm);
                    k = nrejc * (nrejc + 1) / 2 + 1;
                    for (np = 1; np <= nPerm; np++)
                    {
                        XPerm(genomeData, px, rnd);
                        pstat = CBSTStatistic.HTMaxP(hk, totalSumOfSquares, px, sx, al0, isBinary);
                        if (ostat <= pstat)
                        {
                            nrej++;
                            k++;
                        }
                        if (nrej > nrejc) { return; }
                        if (np >= sbdry[k - 1]) { break; }
                    }
                }
                else
                {
                    nrejc = (int)(cutoffPValue * nPerm);
                    k = nrejc * (nrejc + 1) / 2 + 1;
                    for (np = 1; np <= nPerm; np++)
                    {
                        XPerm(genomeData, px, rnd);
                        pstat = CBSTStatistic.TMaxP(totalSumOfSquares, px, sx, al0, isBinary);
                        if (ostat <= pstat)
                        {
                            nrej++;
                            k++;
                        }
                        if (nrej > nrejc) { return; }
                        if (np >= sbdry[k - 1]) { break; }
                    }
                }
            }
            // 200
            if (iseg[1] == genomeData.Length) // The second change point is the right boundary
            {
                nChangePoints = 1;
                iChangePoint[0] = iseg[0];
            }
            else
            {
                if (iseg[0] == 0) // The first change point is the left boundary
                {
                    nChangePoints = 1;
                    iChangePoint[0] = iseg[1];
                }
                else
                {
                    l = 0;
                    n1 = iseg[0];
                    n12 = iseg[1];
                    n2 = n12 - n1;
                    // |-- n1 = iseg[0] --|-- n2 = n12 - n1 --|
                    // |--            n12 = iseg[1]         --|
                    tPValue = CBSTStatistic.TPermP(n1, n2, n12, genomeData, l, px, nPerm, rnd);
                    if (tPValue <= cutoffPValue)
                    {
                        nChangePoints = 1;
                        iChangePoint[0] = iseg[0];
                    }
                    l = iseg[0];
                    n12 = genomeData.Length - iseg[0];
                    n2 = genomeData.Length - iseg[1];
                    n1 = n12 - n2;
                    // |-- n1 = n12 - n2 --|-- n2 = n - iseg[1] --|
                    // |--         n12 = n - iseg[0]            --|
                    tPValue = CBSTStatistic.TPermP(n1, n2, n12, genomeData, l, px, nPerm, rnd);
                    if (tPValue <= cutoffPValue)
                    {
                        nChangePoints++;
                        iChangePoint[nChangePoints - 1] = iseg[1];
                    }
                }
            }
            // 500
        }

        /// <summary>
        /// Permute x and store the results in px
        /// </summary>
        /// <param name="x">the array to be permuted</param>
        /// <param name="px">the permuted array</param>
        private static void XPerm(double[] x, double[] px, Random rnd)
        {
            int i, j;
            double cc;

            for (i = 0; i < x.Length; i++) { px[i] = x[i]; }

            for (i = x.Length - 1; i >= 0; i--)
            {
                cc = rnd.NextDouble();
                j = (int)(cc * (i + 1));
                j = (j > i) ? i : j;
                Helper.Swap<double>(ref px[i], ref px[j]);
            }
        }

        public static double TrimmedVariance(IDictionary<string, double[]> scoresByChr, double trim = 0.025)
        {
            int n = 0;
            foreach (string chr in scoresByChr.Keys)
            {
                n += scoresByChr[chr].Length;
            }
            double[] diff = new double[n - 1];
            int i = 0;
            double last = Double.NaN;
            foreach (string chr in scoresByChr.Keys)
            {
                if (scoresByChr[chr].Length <= 0) { continue; }
                if (i > 0)
                {
                    diff[i] = scoresByChr[chr][0] - last;
                    i++;
                }
                Array.Copy(Helper.Diff(scoresByChr[chr]), 0, diff, i, scoresByChr[chr].Length - 1);
                i += (scoresByChr[chr].Length - 1);
                last = scoresByChr[chr][scoresByChr[chr].Length - 1];
            }
            int nKeep = Convert.ToInt32(Math.Round((1 - 2 * trim) * (n - 1)));
            // R code: inflfact(trim)*sum((sort(abs(diff(genomdat)))[1:n.keep])^2 / (2*n.keep))

            Helper.InplaceAbs(diff);
            Array.Sort(diff);

            return ChangePoint.InflationFactor(trim) * Helper.PartialSumOfPowers(diff, 2, 0, nKeep) / (2 * nKeep);
        }


        /// <summary>
        /// Truncated N(0, 1) to (x1, x2), where P(X <= x1) = trim and P(X <= x2) = 1 - trim
        /// </summary>
        /// <param name="trim">tail probability</param>
        /// <returns>approximately (E[X^2] = 1) / Et[X^2]</returns>
        private static double InflationFactor(double trim)
        {
            var norm = new Normal(); // N(0, 1)
            double a = norm.InverseCumulativeDistribution(1 - trim);
            double step = 2 * a / 10000;
            double[] x1s = Helper.Seq(-a + step / 2, a - step / 2, 10000);
            // Truncated N(0, 1) to P(X <= x) = trim and P(X <= x) = 1 - trim
            double eX2 = 0.0;
            foreach (double x1 in x1s) { eX2 += (x1 * x1) * norm.Density(x1); }
            eX2 = eX2 * step / (1 - 2 * trim);
            // eX2 now approximates Et[X^2]: E[X^2] of the truncated N(0, 1)
            return 1 / eX2; // approx (E[X^2] = 1) / Et[X^2]
            // == 1 / (1 + (-a * dnorm(-a) - a * dnorm(a)) / (1 - 2*trim) - ((dnorm(-1) - dnorm(a)) / (1 - 2* trim))^2)
            // According to http://en.wikipedia.org/wiki/Truncated_normal_distribution
        }

    }
}