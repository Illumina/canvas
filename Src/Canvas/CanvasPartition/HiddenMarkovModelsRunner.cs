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

        public HiddenMarkovModelsRunner(string commonCNVs, int minSize = 10, int nHiddenStates = 5)
        {
            _commonCnVs = commonCNVs;
            _nHiddenStates = nHiddenStates;
            _minSize = minSize;
        }

        public Dictionary<string, Segmentation.Segment[]> Run(Segmentation segmentation) {
            Dictionary<string, List<GenomicBin>> commonCNVintervals = null;
            if (_commonCnVs != null)
            {
                commonCNVintervals = CanvasCommon.Utilities.LoadBedFile(_commonCnVs);
                CanvasCommon.Utilities.SortAndOverlapCheck(commonCNVintervals, _commonCnVs);
            }
            Dictionary<string, Segmentation.Segment[]> segmentByChr = new Dictionary<string, Segmentation.Segment[]>();

            var cts = new CancellationTokenSource();
            Parallel.ForEach(
                segmentation.ScoresByChr.Keys,
                new ParallelOptions
                {
                    CancellationToken = cts.Token,
                    MaxDegreeOfParallelism = Environment.ProcessorCount,
                    TaskScheduler = TaskScheduler.Default
                },
                chr =>
                {
                    List<int> breakpointsForward = new List<int>();
                    List<int> breakpointsReverse = new List<int>();

                    int length = segmentation.ScoresByChr[chr].Count;
                    if (length > _minSize)
                    {
                        List<double> haploidMeans = new List<double>(_nHiddenStates);
                        List<MultivariateNegativeBinomial> negativeBinomialDistributions = InitializeNegativeBinomialEmission(segmentation.ScoresByChr, chr, _nHiddenStates, haploidMeans);
                        HiddenMarkovModel hmm = new HiddenMarkovModel(segmentation.ScoresByChr[chr], negativeBinomialDistributions, haploidMeans);
                        Console.WriteLine($"{DateTime.Now} Launching HMM task for chromosome {chr}");
                        // forward sequence
                        hmm.FindMaximalLikelihood(segmentation.ScoresByChr[chr], chr);
                        List<int> hiddenStatesOriginal = hmm.BestPathViterbi(segmentation.ScoresByChr[chr], segmentation.StartByChr[chr], haploidMeans);
                        //reverse sequence
                        // segmentation.ScoresByChr[chr].Reverse();
                        // hmm.FindMaximalLikelihood(segmentation.ScoresByChr[chr], chr);
                        // List<int> hiddenStatesReverse = hmm.BestPathViterbi(segmentation.ScoresByChr[chr], segmentation.StartByChr[chr], haploidMeans);
                        // segmentation.ScoresByChr[chr].Reverse();

                        // hiddenStatesReverse.Reverse();
                        Console.WriteLine($"{DateTime.Now} Completed HMM task for chromosome {chr}");

                        breakpointsForward.Add(0);
                        // breakpointsReverse.Add(0);
                        for (int i = 1; i < length; i++)
                        {
                            if (hiddenStatesOriginal[i] - hiddenStatesOriginal[i - 1] != 0)
                                breakpointsForward.Add(i);
                            // if (hiddenStatesReverse[i] - hiddenStatesReverse[i - 1] != 0)
                            //  breakpointsReverse.Add(i);
                        }

                        var breakpoints = breakpointsForward; //MergeBreakpoint(breakpointsForward, breakpointsReverse);

                        if (_commonCnVs != null)
                        {
                            if (commonCNVintervals.ContainsKey(chr))
                            {
                                List<GenomicBin> remappedCommonCNVintervals = Segmentation.RemapCommonRegions(commonCNVintervals[chr], segmentation.StartByChr[chr], segmentation.EndByChr[chr]);
                                List<int> oldbreakpoints = breakpoints;
                                breakpoints = Segmentation.OverlapCommonRegions(oldbreakpoints, remappedCommonCNVintervals);
                            }
                        }

                        var segments = segmentation.DeriveSegments(breakpoints, length, chr);

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

        public List<MultivariateGaussianDistribution> InitializeEmission(Dictionary<string, List<List<double>>> data ,string chromosome, int nHiddenStates)
        {
            int nDimensions = data[chromosome].First().Count;
            List<double> haploidMean = new List<double>(nDimensions);
            List<double> standardDeviation = new List<double>(nDimensions);
            List<MultivariateGaussianDistribution> tmpDistributions = new List<MultivariateGaussianDistribution>();

            for (int dimension = 0; dimension < nDimensions; dimension++)
            {
                double meanHolder = 0;
                foreach (List<double> datapoint in data[chromosome])
                    meanHolder += datapoint[dimension];
                haploidMean.Add(meanHolder / data[chromosome].Count /2.0);
                standardDeviation.Add(CanvasCommon.Utilities.StandardDeviation(data[chromosome].Select(x => x[dimension]).ToList()));
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
                MultivariateGaussianDistribution tmpDistribution = new MultivariateGaussianDistribution(tmpMean, tmpSds);
                tmpDistributions.Add(tmpDistribution);
            }

            // remove outliers 
            double maxThreshold = tmpDistributions.Last().Mean().Max() * 1.2;
            RemoveOutliers(data, maxThreshold, nDimensions);

            return tmpDistributions;
        }

        public List<MultivariatePoissonDistribution> InitializePoissonEmission(Dictionary<string, List<List<double>>> data, string chromosome, int nHiddenStates, List<double> haploidMean)
        {
            int nDimensions = data[chromosome].First().Count;
            List<MultivariatePoissonDistribution> tmpDistributions = new List<MultivariatePoissonDistribution>();

            for (int dimension = 0; dimension < nDimensions; dimension++)
            {
                double meanHolder = 0;
                foreach (List<double> datapoint in data[chromosome])
                    meanHolder += datapoint[dimension];
                haploidMean.Add(meanHolder / data[chromosome].Count / 2.0);
            }

            for (int CN = 0; CN < nHiddenStates; CN++)
            {
                double scaler = 0.8;
                Vector<double> tmpMean = Vector<double>.Build.Dense(haploidMean.Select(x => Math.Max(CN, 0.1) * scaler * x).ToArray());
                // if few hidden states, increase the last CN state by diploid rather than haploid increment
                if (nHiddenStates < 5 && CN - 1 == nHiddenStates)
                    tmpMean = Vector<double>.Build.Dense(haploidMean.Select(x => Math.Max(CN, 0.5) * x + x).ToArray());
                MultivariatePoissonDistribution tmpDistribution = new MultivariatePoissonDistribution(tmpMean.ToList());
                tmpDistributions.Add(tmpDistribution);
            }

            // remove outliers 
            double maxThreshold = tmpDistributions.Last().Mean().Max() * 1.2;
            RemoveOutliers(data, maxThreshold, nDimensions);

            return tmpDistributions;
        }


        public List<MultivariateNegativeBinomial> InitializeNegativeBinomialEmission(Dictionary<string, List<List<double>>> data, string chromosome, int nHiddenStates, List<double> haploidMean)
        {
            int nDimensions = data[chromosome].First().Count;         
            var variance = new List<double>(nDimensions);
            List<MultivariateNegativeBinomial> tmpDistributions = new List<MultivariateNegativeBinomial>();

            for (int dimension = 0; dimension < nDimensions; dimension++)
            {
                double meanHolder = 0;
                foreach (List<double> datapoint in data[chromosome])
                    meanHolder += datapoint[dimension];
                haploidMean.Add(meanHolder / data[chromosome].Count / 2.0);
                variance.Add(CanvasCommon.Utilities.Variance(data[chromosome].Select(x => x[dimension]).ToList()));
            }

            for (int CN = 0; CN < nHiddenStates; CN++)
            {
                Vector<double> tmpMean = Vector<double>.Build.Dense(haploidMean.Select(x => Math.Max(CN, 0.1) * x).ToArray());
                // if few hidden states, increase the last CN state by diploid rather than haploid increment
                if (nHiddenStates < 5 && CN - 1 == nHiddenStates)
                    tmpMean = Vector<double>.Build.Dense(haploidMean.Select(x => Math.Max(CN, 0.5) * x + x).ToArray());
                MultivariateNegativeBinomial tmpDistribution = new MultivariateNegativeBinomial(tmpMean.ToList(), variance);
                tmpDistributions.Add(tmpDistribution);
            }

            // remove outliers 
            double maxThreshold = tmpDistributions.Last().Mean().Max() * 1.2;
            RemoveOutliers(data, chromosome, maxThreshold, nDimensions);

            return tmpDistributions;
        }

        private static void RemoveOutliers(Dictionary<string, List<List<double>>> data, string chr, double maxThreshold, int nDimensions)
        {           
                for (int length = 0; length < data[chr].Count; length++)
                    for (int dimension = 0; dimension < nDimensions; dimension++)
                        data[chr][length][dimension] = data[chr][length][dimension] > maxThreshold
                            ? maxThreshold
                            : data[chr][length][dimension];
        }
    }
}