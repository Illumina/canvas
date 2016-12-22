using System;
using System.Collections.Generic;
using System.Threading;
using Isas.Shared.Utilities;
using MathNet.Numerics.Random;

namespace CanvasPartition
{
    class CBSRunner
    {
        public static readonly double DefaultAlpha = 0.01;
        private SegmentSplitUndo _undoMethod = SegmentSplitUndo.None;
        private double _alpha;
        private int _maxInterBinDistInSegment;

        public CBSRunner(int maxInterBinDistInSegment, SegmentSplitUndo undoMethod, double alpha)
        {
            _maxInterBinDistInSegment = maxInterBinDistInSegment;
            _undoMethod = undoMethod;
            _alpha = alpha;
        }


        /// <summary>
        /// CBS: circular binary segmentation porting the R function segment in DNAcopy
        /// </summary>
        /// <param name="alpha">Now in this.Alpha</param>
        /// <param name="nPerm"></param>
        /// <param name="pMethod">"hybrid" or "perm"</param>
        /// <param name="minWidth"></param>
        /// <param name="kMax"></param>
        /// <param name="nMin"></param>
        /// <param name="eta"></param>
        /// <param name="sbdry"></param>
        /// <param name="trim"></param>
        /// <param name="undoSplit">"none" or "prune" or "sdundo"; now in this.UndoMethod</param>
        /// <param name="undoPrune"></param>
        /// <param name="undoSD"></param>
        /// <param name="verbose"></param>
        public Dictionary<string, Segmentation.Segment[]> Run(Segmentation segmentation, uint nPerm = 10000, string pMethod = "hybrid", int minWidth = 2, int kMax = 25,
            uint nMin = 200, double eta = 0.05, uint[] sbdry = null, double trim = 0.025,
            double undoPrune = 0.05, double undoSD = 3, int verbose = 1)
        {
            if (minWidth < 2 || minWidth > 5)
            {
                Console.Error.WriteLine("Minimum segment width should be between 2 and 5");
                Environment.Exit(1);
            }
            if (nMin < 4 * kMax)
            {
                Console.Error.WriteLine("nMin should be >= 4 * kMax");
                Environment.Exit(1);
            }
            if (sbdry == null)
            {
                GetBoundary.ComputeBoundary(nPerm, this._alpha, eta, out sbdry);
            }

            Dictionary<string, int[]> inaByChr = new Dictionary<string, int[]>();
            Dictionary<string, double[]> finiteScoresByChr = new Dictionary<string, double[]>();

            List<ThreadStart> tasks = new List<ThreadStart>();
            foreach (KeyValuePair<string, double[]> scoreByChrKVP in segmentation.ScoreByChr)
            {
                tasks.Add(new ThreadStart(() =>
                {
                    string chr = scoreByChrKVP.Key;
                    int[] ina;
                    Helper.GetFiniteIndices(scoreByChrKVP.Value, out ina); // not NaN, -Inf, Inf

                    double[] scores;
                    if (ina.Length == scoreByChrKVP.Value.Length)
                    {
                        scores = scoreByChrKVP.Value;
                    }
                    else
                    {
                        Helper.ExtractValues<double>(scoreByChrKVP.Value, ina, out scores);
                    }

                    lock (finiteScoresByChr)
                    {
                        finiteScoresByChr[chr] = scores;
                        inaByChr[chr] = ina;
                    }

                }));
            }
            Isas.Shared.Utilities.Utilities.DoWorkParallelThreads(tasks);

            // Quick sanity-check: If we don't have any segments, then return a dummy result.
            int n = 0;
            foreach (var list in finiteScoresByChr.Values)
            {
                n += list.Length;
            }
            if (n == 0)
            {
                return new Dictionary<string, Segmentation.Segment[]>();
            }

            double trimmedSD = Math.Sqrt(ChangePoint.TrimmedVariance(finiteScoresByChr, trim: trim));

            Dictionary<string, Segmentation.Segment[]> segmentByChr = new Dictionary<string, Segmentation.Segment[]>();

            // when parallelizing we need an RNG for each chromosome to get deterministic results
            Random seedGenerator = new MersenneTwister(0);
            Dictionary<string, Random> perChromosomeRandom = new Dictionary<string, Random>();
            foreach (string chr in segmentation.ScoreByChr.Keys)
            {
                perChromosomeRandom[chr] = new MersenneTwister(seedGenerator.NextFullRangeInt32(), true);
            }

            tasks = new List<ThreadStart>();
            foreach (string chr in segmentation.ScoreByChr.Keys)
            {
                tasks.Add(new ThreadStart(() =>
                {
                    int[] ina = inaByChr[chr];
                    int[] lengthSeg;
                    double[] segmentMeans;
                    ChangePoint.ChangePoints(segmentation.ScoreByChr[chr], sbdry, out lengthSeg, out segmentMeans, perChromosomeRandom[chr],
                        dataType: "logratio", alpha: this._alpha, nPerm: nPerm,
                        pMethod: pMethod, minWidth: minWidth, kMax: kMax, nMin: nMin, trimmedSD: trimmedSD,
                        undoSplits: this._undoMethod, undoPrune: undoPrune, undoSD: undoSD, verbose: verbose);

                    Segmentation.Segment[] segments = new Segmentation.Segment[lengthSeg.Length];
                    int cs1 = 0, cs2 = -1; // cumulative sum
                    for (int i = 0; i < lengthSeg.Length; i++)
                    {
                        cs2 += lengthSeg[i];
                        int start = ina[cs1];
                        int end = ina[cs2];
                        segments[i] = new Segmentation.Segment();
                        segments[i].start = segmentation.StartByChr[chr][start]; // Genomic start
                        segments[i].end = segmentation.EndByChr[chr][end]; // Genomic end
                        cs1 += lengthSeg[i];
                    }

                    lock (segmentByChr)
                    {
                        segmentByChr[chr] = segments;
                    }
                }));
            }
            
            Isas.Shared.Utilities.Utilities.DoWorkParallelThreads(tasks);
            // segmentation.SegmentationResults = new Segmentation.GenomeSegmentationResults(segmentByChr);
            Console.WriteLine("{0} Completed CBS tasks", DateTime.Now);
            Console.WriteLine("{0} Segmentation results complete", DateTime.Now);
            return segmentByChr;
        }
    }
}