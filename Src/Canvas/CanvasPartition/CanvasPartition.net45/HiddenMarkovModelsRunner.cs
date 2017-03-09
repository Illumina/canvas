using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using CanvasCommon;
using MathNet.Numerics.LinearAlgebra;

namespace CanvasPartition
{
    class HiddenMarkovModelsRunner
    {
        private int _minSize;
        private int _nHiddenStates;
        private string _commonCnVs;
        private int _nSamples;


        public HiddenMarkovModelsRunner(string commonCNVs, int nSamples, int minSize = 10, int nHiddenStates = 5)
        {
            _commonCnVs = commonCNVs;
            _nHiddenStates = nHiddenStates;
            _minSize = minSize;
            _nSamples = nSamples;
        }

        public Dictionary<string, Segmentation.Segment[]> Run(List<Segmentation> segmentation) {
            Dictionary<string, List<SampleGenomicBin>> commonCNVintervals = null;
            if (_commonCnVs != null)
            {
                commonCNVintervals = CanvasCommon.Utilities.LoadBedFile(_commonCnVs);
                CanvasCommon.Utilities.SortAndOverlapCheck(commonCNVintervals, _commonCnVs);
            }
            Dictionary<string, Segmentation.Segment[]> segmentByChr = new Dictionary<string, Segmentation.Segment[]>();

            var cts = new CancellationTokenSource();
            Parallel.ForEach(
                segmentation.First().ScoreByChr.Keys,
                new ParallelOptions
                {
                    CancellationToken = cts.Token,
                    MaxDegreeOfParallelism = Environment.ProcessorCount,
                    TaskScheduler = TaskScheduler.Default
                },
                chr =>
                {
                    List<int> breakpoints = new List<int>();
                    int length = segmentation.First().ScoreByChr[chr].Length;
                    var startByChr = segmentation.First().StartByChr[chr];
                    var endByChr = segmentation.First().EndByChr[chr];
                    List<List<double>> multiSampleCoverage = new List<List<double>>(length);
                    for (int i = 0; i < length; i++)
                        multiSampleCoverage.Add(segmentation.Select(x=>x.ScoreByChr[chr][i]).ToList());

                    if (length > _minSize)
                    {
                        List<double> haploidMeans = new List<double>(_nHiddenStates);
                        List<MultivariateNegativeBinomial> negativeBinomialDistributions = InitializeNegativeBinomialEmission(multiSampleCoverage, _nHiddenStates, haploidMeans);
                        HiddenMarkovModel hmm = new HiddenMarkovModel(multiSampleCoverage, negativeBinomialDistributions, haploidMeans);
                        Console.WriteLine($"{DateTime.Now} Launching HMM task for chromosome {chr}");
                        if (_nSamples <= 3)
                        hmm.FindMaximalLikelihood(multiSampleCoverage);
                        List<int> bestPathViterbi = hmm.BestPathViterbi(multiSampleCoverage, startByChr, haploidMeans);
                        Console.WriteLine($"{DateTime.Now} Completed HMM task for chromosome {chr}");

                        breakpoints.Add(0);
                        for (int i = 1; i < length; i++)
                        {
                            if (bestPathViterbi[i] - bestPathViterbi[i - 1] != 0)
                            {
                                breakpoints.Add(i);
                            }
                        }


                        if (_commonCnVs != null)
                        {
                            if (commonCNVintervals.ContainsKey(chr))
                            {
                                List<SampleGenomicBin> remappedCommonCNVintervals = Segmentation.RemapCommonRegions(commonCNVintervals[chr], startByChr, endByChr);
                                List<int> oldbreakpoints = breakpoints;
                                breakpoints = Segmentation.OverlapCommonRegions(oldbreakpoints, remappedCommonCNVintervals);
                            }
                        }

                        var segments = Segmentation.DeriveSegments(breakpoints, length, startByChr, endByChr);

                        lock (segmentByChr)
                        {
                            segmentByChr[chr] = segments;
                        }
                    }
                });

            Console.WriteLine("{0} Completed HMM tasks", DateTime.Now);
            Console.WriteLine("{0} Segmentation results complete", DateTime.Now);
            return segmentByChr;
        }

        private static List<int> MergeBreakpoint(List<int> breakpointsForward, List<int> breakpointsReverse)
        {
            var expandedReverseBreakpoints = breakpointsReverse.Select(x => x - 1).Union(breakpointsReverse.Select(x => x + 1)).Union(breakpointsReverse);
            var truncatedBreakpointsForward = breakpointsForward.Except(expandedReverseBreakpoints);
            var breakpointsMerged = truncatedBreakpointsForward.Union(breakpointsReverse).ToList();
            breakpointsMerged.Sort();
            var breakpointsMergeFiltered = new List<int>(breakpointsMerged.Count);
            breakpointsMergeFiltered.Add(breakpointsMerged.First());
            for (int i = 1; i < breakpointsMerged.Count; i++)
            {
                if (breakpointsMerged[i]!= breakpointsMerged[i-1]+1)
                    breakpointsMergeFiltered.Add(breakpointsMerged[i]);
            }
            return breakpointsMergeFiltered;
        }

        public List<MultivariateGaussianDistribution> InitializeEmission(List<List<double>> data ,string chromosome, int nHiddenStates)
        {
            int nDimensions = _nSamples;
            List<double> haploidMean = new List<double>(nDimensions);
            List<double> standardDeviation = new List<double>(nDimensions);
            var tmpDistributions = new List<MultivariateGaussianDistribution>();

            for (int dimension = 0; dimension < nDimensions; dimension++)
            {
                double meanHolder = data.Sum(datapoint => datapoint[dimension]);
                haploidMean.Add(meanHolder / data.Count /2.0);
                standardDeviation.Add(CanvasCommon.Utilities.StandardDeviation(data.Select(x => x[dimension]).ToList()));
            }

            for (int CN = 0; CN < nHiddenStates; CN++)
            {
                Matrix<double> tmpSds = Matrix<double>.Build.Dense(nDimensions, nDimensions, 0);

                for (int dimension = 0; dimension < nDimensions; dimension++)
                {
                    tmpSds[dimension, dimension] = standardDeviation[dimension];
                }
                Vector<double> tmpMean = Vector<double>.Build.Dense(haploidMean.Select(x => Math.Max(CN, 0.05)*x).ToArray());
                // if few hidden states, increase the last CN state by diploid rather than haploid increment
                if (nHiddenStates < 5 && CN-1 == nHiddenStates)
                    tmpMean = Vector<double>.Build.Dense(haploidMean.Select(x => Math.Max(CN, 0.5) * x + x).ToArray());
                var tmpDistribution = new MultivariateGaussianDistribution(tmpMean, tmpSds);
                tmpDistributions.Add(tmpDistribution);
            }

            // remove outliers 
            double maxThreshold = tmpDistributions.Last().Mean().Max() * 1.2;
            RemoveOutliers(data, maxThreshold);

            return tmpDistributions;
        }

        public List<MultivariatePoissonDistribution> InitializePoissonEmission(List<List<double>> data,  int nHiddenStates)
        {
            int nDimensions = _nSamples;
            List<double> haploidMean = new List<double>(nDimensions);
            var tmpDistributions = new List<MultivariatePoissonDistribution>();

            for (int dimension = 0; dimension < nDimensions; dimension++)
            {
                double meanHolder = data.Sum(datapoint => datapoint[dimension]);
                haploidMean.Add(meanHolder / data.Count / 2.0);
            }

            for (int CN = 0; CN < nHiddenStates; CN++)
            {
                double scaler = 0.8;
                Vector<double> tmpMean = Vector<double>.Build.Dense(haploidMean.Select(x => Math.Max(CN, 0.1) * scaler * x).ToArray());
                // if few hidden states, increase the last CN state by diploid rather than haploid increment
                if (nHiddenStates < 5 && CN - 1 == nHiddenStates)
                    tmpMean = Vector<double>.Build.Dense(haploidMean.Select(x => Math.Max(CN, 0.5) * x + x).ToArray());
                var tmpDistribution = new MultivariatePoissonDistribution(tmpMean.ToList());
                tmpDistributions.Add(tmpDistribution);
            }

            // remove outliers 
            double maxThreshold = tmpDistributions.Last().Mean().Max() * 1.2;
            RemoveOutliers(data, maxThreshold);

            return tmpDistributions;
        }


        public List<MultivariateNegativeBinomial> InitializeNegativeBinomialEmission(List<List<double>> data, int nHiddenStates, List<double> haploidMean)
        {
            int nDimensions = _nSamples;         
            var variance = new List<double>(nDimensions);
            List<MultivariateNegativeBinomial> tmpDistributions = new List<MultivariateNegativeBinomial>();

            for (int dimension = 0; dimension < nDimensions; dimension++)
            {
                double meanHolder = data.Sum(datapoint => datapoint[dimension]);
                haploidMean.Add(meanHolder / data.Count / 2.0);
                variance.Add(CanvasCommon.Utilities.Variance(data.Select(x => x[dimension]).ToList()));
            }

            // remove outliers 
            double maxThreshold = haploidMean.Max() * nHiddenStates;
            RemoveOutliers(data, maxThreshold);
            var maxValues = data.Select(x => Convert.ToInt32(x.Max())).ToList().Max();

            for (int CN = 0; CN < nHiddenStates; CN++)
            {
                Vector<double> tmpMean = Vector<double>.Build.Dense(haploidMean.Select(x => Math.Max(CN, 0.1) * x).ToArray());
                // if few hidden states, increase the last CN state by diploid rather than haploid increment
                if (nHiddenStates < 5 && CN - 1 == nHiddenStates)
                    tmpMean = Vector<double>.Build.Dense(haploidMean.Select(x => Math.Max(CN, 0.5) * x + x).ToArray());
                MultivariateNegativeBinomial tmpDistribution = new MultivariateNegativeBinomial(tmpMean.ToList(), variance, maxValues + 10);
                tmpDistributions.Add(tmpDistribution);
            }

            return tmpDistributions;
        }

        private static void RemoveOutliers(List<List<double>> data, double maxThreshold)
        {
            int nDimensions = data.First().Count;
            foreach (List<double> sample in data)
                for (int dimension = 0; dimension < nDimensions; dimension++)
                    sample[dimension] = sample[dimension] > maxThreshold
                        ? maxThreshold
                        : sample[dimension];
        }
    }
}