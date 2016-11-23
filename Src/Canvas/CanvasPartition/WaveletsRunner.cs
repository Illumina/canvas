using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using CanvasCommon;

namespace CanvasPartition
{
    class WaveletsRunner
    {
        private WaveletsRunnerParams _parameters;
        public static readonly double DefaultMadFactor = 2.0; // default MAD factor for Wavelets

        public class WaveletsRunnerParams
        {
            public bool IsGermline { get; }
            public string CommonCnVs { get; }
            public double ThresholdLower { get; }
            public double ThresholdUpper { get; }
            public double MadFactor { get; }
            public int MinSize { get; }
            public int Verbose { get; }
            public bool IsSmallPedegree { get; }

            public WaveletsRunnerParams(bool isGermline, string commonCNVs, double thresholdLower = 5, double thresholdUpper = 80, double madFactor = 2, int minSize = 10, int verbose = 1, bool isSmallPedegree = false)
            {
                IsGermline = isGermline;
                CommonCnVs = commonCNVs;
                ThresholdLower = thresholdLower;
                ThresholdUpper = thresholdUpper;
                MadFactor = madFactor;
                MinSize = minSize;
                Verbose = verbose;
                IsSmallPedegree = isSmallPedegree;
            }
        }

        public WaveletsRunner(WaveletsRunnerParams waveletsRunnerParams)
        {
            _parameters = waveletsRunnerParams;
        }

        /// <summary>
        /// Wavelets: unbalanced HAAR wavelets segmentation 
        /// </summary>
        /// <param name="threshold">wavelets coefficient threshold</param>
        public Dictionary<string, Segmentation.Segment[]> Run(Segmentation segmentation)
        {
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
            int n = finiteScoresByChr.Values.Sum(list => list.Length);

            if (n == 0)
                return new Dictionary<string, Segmentation.Segment[]>();           

            Dictionary<string, Segmentation.Segment[]> segmentByChr = new Dictionary<string, Segmentation.Segment[]>();

            // load common CNV segments
            Dictionary<string, List<SampleGenomicBin>> commonCNVintervals = null;
            if (_parameters.CommonCnVs != null)
            {
                commonCNVintervals = CanvasCommon.Utilities.LoadBedFile(_parameters.CommonCnVs);
                CanvasCommon.Utilities.SortAndOverlapCheck(commonCNVintervals, _parameters.CommonCnVs);
            }

            tasks = new List<ThreadStart>();
            foreach (string chr in segmentation.ScoreByChr.Keys)
            {
                tasks.Add(new ThreadStart(() =>
                {
                    List<int> breakpoints = new List<int>();
                    int sizeScoreByChr = segmentation.ScoreByChr[chr].Length;
                    if (sizeScoreByChr > _parameters.MinSize)
                    {
                        WaveletSegmentation.HaarWavelets(segmentation.ScoreByChr[chr], _parameters.ThresholdLower, _parameters.ThresholdUpper,
                            breakpoints, _parameters.IsGermline, madFactor: _parameters.MadFactor);
                    }

                    if (_parameters.CommonCnVs != null)
                    {
                        if (commonCNVintervals.ContainsKey(chr))
                        {
                            List <SampleGenomicBin> remappedCommonCNVintervals = Segmentation.RemapCommonRegions(commonCNVintervals[chr], segmentation.StartByChr[chr], segmentation.EndByChr[chr]);
                            List <int> oldbreakpoints = breakpoints;
                            breakpoints = Segmentation.OverlapCommonRegions(oldbreakpoints, remappedCommonCNVintervals);
                        }
                    }

                    var segments = Segmentation.DeriveSegments(breakpoints, sizeScoreByChr, segmentation.StartByChr[chr], segmentation.EndByChr[chr]);

                    lock (segmentByChr)
                    {
                        segmentByChr[chr] = segments;
                    }
                }));

            }
            Console.WriteLine("{0} Launching wavelet tasks", DateTime.Now);
            Isas.Shared.Utilities.Utilities.DoWorkParallelThreads(tasks);
            Console.WriteLine("{0} Completed wavelet tasks", DateTime.Now);
            Console.WriteLine("{0} Segmentation results complete", DateTime.Now);
            return segmentByChr;
        }
    }
}