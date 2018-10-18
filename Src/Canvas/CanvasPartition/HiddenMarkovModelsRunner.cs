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
        private readonly int _nSamples;

        public HiddenMarkovModelsRunner(int nSamples, int minSize = 10, int nHiddenStates = 5)
        {
            _nHiddenStates = nHiddenStates;
            _minSize = minSize;
            _nSamples = nSamples;
        }

        public Dictionary<string, SegmentationInput.Segment[]> Run(List<SegmentationInput> segmentation, bool isPerSample)
        {
            var segmentByChr = new Dictionary<string, SegmentationInput.Segment[]>();

            var cts = new CancellationTokenSource();

            // Compute whole-genome median and inter-quartile-range-based pseudo-variance for each sample;
            // it would be better to exclude regions that are not diploid, and we should really be
            // using a different variance for each copy number, but using these values is better than
            // using the per-chromosome mean and variance, which have the following problems:
            // - chromosomes with a lot of outliers can get a very high variance
            // - chromosomes that have a whole-chromosome CNV or a CNV that affects a lot of the chromosome
            //   can have problematic estimates
            var medians = new List<double>();
            var pseudoVariances = new List<double>();
            foreach (var singleSampleSegmentation in segmentation)
            {
                var cvgVals = new List<float>();
                foreach (var chr in singleSampleSegmentation.CoverageInfo.CoverageByChr.Keys)
                {
                    cvgVals.AddRange(singleSampleSegmentation.CoverageInfo.CoverageByChr[chr].Select(x => (float)x));
                }
                var quartiles = CanvasCommon.Utilities.Quartiles(cvgVals);
                medians.Add(quartiles.Item2);
                var iqr = quartiles.Item3 - quartiles.Item1;
                pseudoVariances.Add(iqr * iqr);
                //Console.WriteLine($"Global estimation of median and pseudovariance: {quartiles.Item2} {iqr * iqr}");
            }
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
                        multiSampleCoverage.Add(segmentation.Select(x => x.CoverageInfo.CoverageByChr[chr][i]).ToList());

                    if (length > _minSize)
                    {

                        var haploidMeans = new List<double>(_nHiddenStates);
                        var negativeBinomialDistributions = isPerSample ?
                            InitializeNegativeBinomialEmission(multiSampleCoverage, _nHiddenStates, haploidMeans, medians, pseudoVariances)
                          : InitializeNegativeBinomialEmission(multiSampleCoverage, _nHiddenStates, haploidMeans, null, null);
                        //for (int j = 0; j < 1; j++)
                        //    for (int i = 0; i < 190; i++)
                        //    {
                        //        Console.WriteLine($"NegBin smp {j} count {i}: {negativeBinomialDistributions[0].Probability(j, i)} {negativeBinomialDistributions[1].Probability(j, i)} {negativeBinomialDistributions[2].Probability(j, i)} {negativeBinomialDistributions[3].Probability(j, i)} {negativeBinomialDistributions[4].Probability(j, i)}");
                        //    }
                        var hmm = new HiddenMarkovModel(multiSampleCoverage, negativeBinomialDistributions, haploidMeans, isPerSample);
                        Console.WriteLine($"{DateTime.Now} Launching HMM task for chromosome {chr}");
                        //if (_nSamples == 1)
                        //    hmm.FindMaximalLikelihood(multiSampleCoverage);
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

        public List<MultivariateNegativeBinomial> InitializeNegativeBinomialEmission(List<List<double>> data, int nHiddenStates, List<double> haploidMean, List<double> medians = null, List<double> pseudoVariances = null)
        {
            int nDimensions = _nSamples;
            var variance = new List<double>(nDimensions);
            var tmpDistributions = new List<MultivariateNegativeBinomial>();

            for (int dimension = 0; dimension < nDimensions; dimension++)
            {
                double median = Math.Max(1d, CanvasCommon.Utilities.Median(data.Select(datapoint => datapoint[dimension])));
                if (medians == null)
                {
                    haploidMean.Add(median / 2.0);
                    variance.Add(CanvasCommon.Utilities.Variance(data.Select(x => x[dimension]).ToList()));
                }
                else
                {
                    haploidMean.Add(medians[dimension] / 2.0);
                    variance.Add(pseudoVariances[dimension]);
                }
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
                //var tmpVar = pseudoVariances == null ?
                //    variance
                //  : variance.Select(x => (Math.Max(CN, 0.5) / 2d) * (Math.Max(CN, 0.5) / 2d) * x);
                var tmpVar = variance;
                var tmpDistribution = new MultivariateNegativeBinomial(tmpMean.ToList(), tmpVar.ToList(), maxValues + 10);
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
