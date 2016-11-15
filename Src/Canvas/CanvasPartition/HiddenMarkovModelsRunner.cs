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

        public HiddenMarkovModelsRunner(string commonCNVs, int minSize = 10, int nHiddenStates = 6)
        {
            _commonCnVs = commonCNVs;
            _nHiddenStates = nHiddenStates;
            _minSize = minSize;
        }

        public Dictionary<string, Segmentation.Segment[]> Run(List<Segmentation> segmentation) {
            Dictionary<string, List<SampleGenomicBin>> commonCNVintervals = null;
            if (_commonCnVs != null)
            {
                commonCNVintervals = Utilities.LoadBedFile(_commonCnVs);
                Utilities.SortAndOverlapCheck(commonCNVintervals, _commonCnVs);
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
                    List<List<double>> multiSampleCoverage = segmentation.Select(x => x.ScoreByChr[chr].ToList()).ToList();
                    if (length > _minSize)
                    {
                        List<MultivariatePoissonDistribution> gaussianDistribution = InitializePoissonEmission(multiSampleCoverage, _nHiddenStates);
                        HiddenMarkovModel hmm = new HiddenMarkovModel(multiSampleCoverage, gaussianDistribution);
                        Console.WriteLine($"{DateTime.Now} Launching HMM task for chromosome {chr}");
                        hmm.FindMaximalLikelihood(multiSampleCoverage);
                        List<int> hiddenStates = hmm.BestPathViterbi(multiSampleCoverage);
                        Console.WriteLine($"{DateTime.Now} Completed HMM task for chromosome {chr}");

                        breakpoints.Add(0);
                        for (int i = 1; i < length; i++)
                        {
                            if (hiddenStates[i] - hiddenStates[i - 1] != 0)
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

        public List<MultivariateGaussianDistribution> InitializeEmission(List<List<double>> data, int nHiddenStates)
        {
            int nDimensions = data.First().Count;
            List<double> haploidMean = new List<double>(nDimensions);
            List<double> standardDeviation = new List<double>(nDimensions);
            List<MultivariateGaussianDistribution> tmpDistributions = new List<MultivariateGaussianDistribution>();

            for (int dimension = 0; dimension < nDimensions; dimension++)
            {
                double meanHolder = 0;
                foreach (List<double> datapoint in data)
                    meanHolder += datapoint[dimension];
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
                MultivariateGaussianDistribution tmpDistribution = new MultivariateGaussianDistribution(tmpMean, tmpSds);
                tmpDistributions.Add(tmpDistribution);
            }

            // remove outliers 
            double maxThreshold = tmpDistributions.Last().Mean().Max() * 1.2;
            return tmpDistributions;
        }

        public List<MultivariatePoissonDistribution> InitializePoissonEmission(List<List<double>> data,  int nHiddenStates)
        {
            int nDimensions = data.First().Count;
            List<double> haploidMean = new List<double>(nDimensions);
            List<MultivariatePoissonDistribution> tmpDistributions = new List<MultivariatePoissonDistribution>();

            for (int dimension = 0; dimension < nDimensions; dimension++)
            {
                double meanHolder = 0;
                foreach (List<double> datapoint in data)
                    meanHolder += datapoint[dimension];
                haploidMean.Add(meanHolder / data.Count / 2.0);
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

        private static void RemoveOutliers(List<List<double>> data, double maxThreshold, int nDimensions)
        {           
                for (int length = 0; length < data.Count; length++)
                    for (int dimension = 0; dimension < nDimensions; dimension++)
                        data[length][dimension] = data[length][dimension] > maxThreshold
                            ? maxThreshold
                            : data[length][dimension];
        }
    }
}