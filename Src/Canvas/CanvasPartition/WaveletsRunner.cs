using System;
using System.Collections.Generic;
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
            private bool _isGermline;
            private string _commonCnVs;
            private double _thresholdLower;
            private double _thresholdUpper;
            private double _madFactor;
            private int _minSize;
            private int _verbose;
            private bool _isSmallPedegree;

            public WaveletsRunnerParams(bool isGermline, string commonCNVs, double thresholdLower = 5, double thresholdUpper = 80, double madFactor = 2, int minSize = 10, int verbose = 1, bool isSmallPedegree = false)
            {
                _isGermline = isGermline;
                _commonCnVs = commonCNVs;
                _thresholdLower = thresholdLower;
                _thresholdUpper = thresholdUpper;
                _madFactor = madFactor;
                _minSize = minSize;
                _verbose = verbose;
                _isSmallPedegree = isSmallPedegree;
            }


            public bool IsGermline
            {
                get { return _isGermline; }
            }

            public string CommonCnVs
            {
                get { return _commonCnVs; }
            }

            public double ThresholdLower
            {
                get { return _thresholdLower; }
            }

            public double ThresholdUpper
            {
                get { return _thresholdUpper; }
            }

            public double MadFactor
            {
                get { return _madFactor; }
            }

            public int MinSize
            {
                get { return _minSize; }
            }

            public int Verbose
            {
                get { return _verbose; }
            }

            public bool IsSmallPedegree
            {
                get { return _isSmallPedegree; }
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
            int n = 0;
            foreach (var list in finiteScoresByChr.Values)
            {
                n += list.Length;
            }
            if (n == 0)
            {
                return new Dictionary<string, Segmentation.Segment[]>();
            }

            Dictionary<string, Segmentation.Segment[]> segmentByChr = new Dictionary<string, Segmentation.Segment[]>();


            // load common CNV segments
            Dictionary<string, List<GenomicBin>> commonCNVintervals = null;
            if (_parameters.CommonCnVs != null)
            {
                commonCNVintervals = Utilities.LoadBedFile(_parameters.CommonCnVs);
                Utilities.SortAndOverlapCheck(commonCNVintervals, _parameters.CommonCnVs);
            }

            tasks = new List<ThreadStart>();
            foreach (string chr in segmentation.ScoreByChr.Keys)
            {
                tasks.Add(new ThreadStart(() =>
                {
                    int[] ina = inaByChr[chr];
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
                            List <GenomicBin> remappedCommonCNVintervals = Segmentation.RemapCommonRegions(commonCNVintervals[chr], segmentation.StartByChr[chr], segmentation.EndByChr[chr]);
                            List <int> oldbreakpoints = breakpoints;
                            breakpoints = Segmentation.OverlapCommonRegions(oldbreakpoints, remappedCommonCNVintervals);
                        }
                    }

                    var segments = segmentation.DeriveSegments(breakpoints, sizeScoreByChr, chr);

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