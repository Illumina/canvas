using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using CanvasCommon;
using Illumina.Common;

namespace CanvasPartition
{
    class WaveletsRunner
    {
        private WaveletsRunnerParams _parameters;
        public static readonly double DefaultMadFactor = 5.0; // default MAD factor for Wavelets

        public class WaveletsRunnerParams
        {
            public bool IsGermline { get; }
            public double EvennessScoreThreshold { get; }
            public string CommonCNVs { get; }
            public double ThresholdLower { get; }
            public double ThresholdUpper { get; }
            public double MadFactor { get; }
            public int MinSize { get; }
            public int Verbose { get; }
            public bool IsSmallPedegree { get; }

            public WaveletsRunnerParams(bool isGermline, string commonCNVs = null, double evennessScoreThreshold = 94.50,
                double thresholdLower = 5, double thresholdLowerMaf = 0.05, double thresholdUpper = 80, double madFactor = 5, int minSize = 10,
                int verbose = 1, bool isSmallPedegree = false)
            {
                IsGermline = isGermline;
                EvennessScoreThreshold = evennessScoreThreshold;
                CommonCNVs = commonCNVs;
                ThresholdLower = thresholdLowerMaf;
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
        public Dictionary<string, SegmentationInput.Segment[]> Run(SegmentationInput segmentationInput, int windowSize)
        {
            double? coverageCV = segmentationInput.GetCoverageVariability(windowSize);
            var factorOfThreeCMADs = segmentationInput.FactorOfThreeCoverageVariabilities(); ;
            try
            {

                double evennessScore = segmentationInput.GetEvennessScore(windowSize);
                if (!segmentationInput.EvennessMetricFile.IsNullOrEmpty())
                    CanvasIO.WriteEvennessMetricToTextFile(segmentationInput.EvennessMetricFile, evennessScore);
            }
            catch (Exception)
            {
                Console.Error.WriteLine("Unable to calculate an evenness score, using coverage for segmentation");
            }

            Dictionary<string, List<int>> adjustedBreakpoints;

            var breakpoints = LaunchWavelets(segmentationInput.CoverageInfo.CoverageByChr, segmentationInput.CoverageInfo.StartByChr,
                segmentationInput.CoverageInfo.EndByChr, coverageCV, factorOfThreeCMADs);
            adjustedBreakpoints = AdjustBreakpoints(segmentationInput.CoverageInfo.CoverageByChr, breakpoints, vafContainingBinsByChr: null);

            var segments = new Dictionary<string, SegmentationInput.Segment[]>();
            foreach (string chr in segmentationInput.VafByChr.Keys)
            {
                segments[chr] = SegmentationInput.DeriveSegments(adjustedBreakpoints[chr], segmentationInput.CoverageInfo.CoverageByChr[chr].Length,
                    segmentationInput.CoverageInfo.StartByChr[chr], segmentationInput.CoverageInfo.EndByChr[chr]);
            }
            return segments;
        }

        public Dictionary<string, List<int>> LaunchWavelets(Dictionary<string, double[]> coverageByChr, Dictionary<string, uint[]> startByChr,
            Dictionary<string, uint[]> endByChr, double? CV, List<double> factorOfThreeCMADs)
        {
            var inaByChr = new Dictionary<string, int[]>();
            var finiteScoresByChr = new Dictionary<string, double[]>();

            var tasks = coverageByChr.Select(scoreByChrKVP => new ThreadStart(() =>
            {
                string chr = scoreByChrKVP.Key;
                Helper.GetFiniteIndices(scoreByChrKVP.Value, out int[] ina); // not NaN, -Inf, Inf

                double[] scores;
                if (ina.Length == scoreByChrKVP.Value.Length)
                    scores = scoreByChrKVP.Value;
                else
                    Helper.ExtractValues<double>(scoreByChrKVP.Value, ina, out scores);

                lock (finiteScoresByChr)
                {
                    finiteScoresByChr[chr] = scores;
                    inaByChr[chr] = ina;
                }
            })).ToList();


            Parallel.ForEach(tasks, task => task.Invoke());
            // Quick sanity-check: If we don't have any segments, then return a dummy result.
            int n = finiteScoresByChr.Values.Sum(list => list.Length);
            if (n == 0)
                return new Dictionary<string, List<int>>();

            var breakpointsByChr = new Dictionary<string, List<int>>();
            tasks = coverageByChr.Keys.Select(chr => new ThreadStart(() =>
            {
                var breakpoints = new List<int>();
                // to cover cases of no SNVs present (i.e. chrY) => chromosome becomes one segment
                int segmentLengthByChr = Math.Max(coverageByChr[chr].Length, 1);
                if (segmentLengthByChr > _parameters.MinSize)
                {
                    WaveletSegmentation.HaarWavelets(coverageByChr[chr], _parameters.ThresholdLower,
                        _parameters.ThresholdUpper,
                        breakpoints, _parameters.IsGermline, _parameters.MadFactor,
                        CV, factorOfThreeCMADs, chr);
                }

                lock (breakpointsByChr)
                {
                    breakpointsByChr[chr] = breakpoints;
                }
            })).ToList();

            Console.WriteLine("{0} Launching wavelet tasks", DateTime.Now);
            Parallel.ForEach(tasks, task => task.Invoke());
            Console.WriteLine("{0} Completed wavelet tasks", DateTime.Now);
            Console.WriteLine("{0} Segmentation results complete", DateTime.Now);
            return breakpointsByChr;
        }

        private Dictionary<string, List<int>> AdjustBreakpoints(Dictionary<string, double[]> binsByChr,
            Dictionary<string, List<int>> breakpoints, Dictionary<string, int[]> vafContainingBinsByChr)
        {
            var adjustedBreakpoints = new Dictionary<string, List<int>>(breakpoints);

            foreach (string chr in binsByChr.Keys)
            {
                if (vafContainingBinsByChr?[chr] != null && vafContainingBinsByChr[chr].Length > 0)
                {
                    if (breakpoints[chr].Max() > vafContainingBinsByChr[chr].Length)
                        throw new Exception(
                            $"breakpoint {breakpoints.Max()} is larger then the zero-based size of the " +
                            $"vafToCoverageIndex '{vafContainingBinsByChr.Count - 1}'");
                    adjustedBreakpoints[chr] = breakpoints[chr].Select(breakpoint => vafContainingBinsByChr[chr][breakpoint]).ToList();
                }
            }
            return adjustedBreakpoints;
        }

        public double[] WaveletMeanSmoother(double[] canvasBins)
        {
            var canvasBinsCopy = new double[canvasBins.Length];
            const int halfWindow = 1;
            for (var index = 0; index < canvasBins.Length; index++)
                canvasBinsCopy[index] = canvasBins[index];
            const double w1 = 1.0 / 3.0;
            const double w2 = 1.0 / 3.0;
            const double w3 = 1.0 / 3.0;

            for (int index = halfWindow; index < canvasBins.Length - halfWindow; index++)
                canvasBinsCopy[index] = canvasBins[index - halfWindow] * w1 + canvasBins[index] * w2 + canvasBins[index + halfWindow] * w3;

            return canvasBinsCopy;
        }
    }
}
