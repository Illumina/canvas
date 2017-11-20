using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;

namespace CanvasPartition
{
    class HiddenMarkovModelsRunner
    {
        private readonly int _minSize;
        private readonly int _nHiddenStates;
        private readonly string _commonCnVs;
        private readonly int _nSamples;

        public HiddenMarkovModelsRunner(string commonCNVs, int nSamples, int minSize = 10, int nHiddenStates = 5)
        {
            _commonCnVs = commonCNVs;
            _nHiddenStates = nHiddenStates;
            _minSize = minSize;
            _nSamples = nSamples;
        }

        public Dictionary<string, SegmentationInput.Segment[]> Run(List<SegmentationInput> segmentation) {
            var segmentByChr = new Dictionary<string, SegmentationInput.Segment[]>();

            var cts = new CancellationTokenSource();
            Parallel.ForEach(
                segmentation.First().CoverageInfo.CoverageByChr.Keys,
                new ParallelOptions
                {
                    CancellationToken = cts.Token,
                    MaxDegreeOfParallelism = Environment.ProcessorCount,
                    TaskScheduler = TaskScheduler.Default
                },
                chr =>
                {
                    var breakpoints = new List<int>();
                    int length = segmentation.First().CoverageInfo.CoverageByChr[chr].Length;
                    var startByChr = segmentation.First().CoverageInfo.StartByChr[chr];
                    var endByChr = segmentation.First().CoverageInfo.EndByChr[chr];
                    var multiSampleCoverage = new List<List<double>>(length);
                    for (int i = 0; i < length; i++)
                        multiSampleCoverage.Add(segmentation.Select(x=>x.CoverageInfo.CoverageByChr[chr][i]).ToList());

                    if (length > _minSize)
                    {
                        var haploidMeans = new List<double>(_nHiddenStates);
                        var negativeBinomialDistributions = InitializeNegativeBinomialEmission(multiSampleCoverage, _nHiddenStates, haploidMeans);
                        var hmm = new HiddenMarkovModel(multiSampleCoverage, negativeBinomialDistributions, haploidMeans);
                        Console.WriteLine($"{DateTime.Now} Launching HMM task for chromosome {chr}");
                        if (_nSamples <= 3)
                            hmm.FindMaximalLikelihood(multiSampleCoverage);
                        var bestPathViterbi = hmm.BestPathViterbi(multiSampleCoverage, startByChr, haploidMeans);
                        Console.WriteLine($"{DateTime.Now} Completed HMM task for chromosome {chr}");

                        breakpoints.Add(0);
                        for (int i = 1; i < length; i++)
                        {
                            if (bestPathViterbi[i] - bestPathViterbi[i - 1] != 0)
                            {
                                breakpoints.Add(i);
                            }
                        }

                        var segments = SegmentationInput.DeriveSegments(breakpoints, length, startByChr, endByChr);

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

        public List<MultivariateGaussianDistribution> InitializeEmission(List<List<double>> data ,string chromosome, int nHiddenStates)
        {
            int nDimensions = _nSamples;
            var haploidMean = new List<double>(nDimensions);
            var standardDeviation = new List<double>(nDimensions);
            var tmpDistributions = new List<MultivariateGaussianDistribution>();

            for (int dimension = 0; dimension < nDimensions; dimension++)
            {
                double meanHolder = data.Sum(datapoint => datapoint[dimension]);
                haploidMean.Add(meanHolder / data.Count /2.0);
                standardDeviation.Add(CanvasCommon.Utilities.StandardDeviation(data.Select(x => x[dimension]).ToList()));
            }

            for (int CN = 0; CN < nHiddenStates; CN++)
            {
                var tmpSds = Matrix<double>.Build.Dense(nDimensions, nDimensions, 0);

                for (int dimension = 0; dimension < nDimensions; dimension++)
                {
                    tmpSds[dimension, dimension] = standardDeviation[dimension];
                }
                var tmpMean = Vector<double>.Build.Dense(haploidMean.Select(x => Math.Max(CN, 0.05)*x).ToArray());
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
            var haploidMean = new List<double>(nDimensions);
            var tmpDistributions = new List<MultivariatePoissonDistribution>();

            for (int dimension = 0; dimension < nDimensions; dimension++)
            {
                double meanHolder = data.Sum(datapoint => datapoint[dimension]);
                haploidMean.Add(meanHolder / data.Count / 2.0);
            }

            for (int CN = 0; CN < nHiddenStates; CN++)
            {
                double scaler = 0.8;
                var tmpMean = Vector<double>.Build.Dense(haploidMean.Select(x => Math.Max(CN, 0.1) * scaler * x).ToArray());
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
            var tmpDistributions = new List<MultivariateNegativeBinomial>();

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
                var tmpMean = Vector<double>.Build.Dense(haploidMean.Select(x => Math.Max(CN, 0.1) * x).ToArray());
                // if few hidden states, increase the last CN state by diploid rather than haploid increment
                if (nHiddenStates < 5 && CN - 1 == nHiddenStates)
                    tmpMean = Vector<double>.Build.Dense(haploidMean.Select(x => Math.Max(CN, 0.5) * x + x).ToArray());
                var tmpDistribution = new MultivariateNegativeBinomial(tmpMean.ToList(), variance, maxValues + 10);
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
