using System;
using System.Collections.Generic;
using System.Linq;
using CanvasCommon;
using Isas.Framework.DataTypes.Maps;

namespace CanvasPedigreeCaller
{
    class CopyNumberLikelihoodCalculator
    {
        private readonly int _maximumCopyNumber;

        public CopyNumberLikelihoodCalculator(int maximumCopyNumber)
        {
            _maximumCopyNumber = maximumCopyNumber;
        }

        /// <summary>
        /// Calculates maximal likelihood for segments without SNV allele ratios. Updated CanvasSegment CopyNumber only. 
        /// Use likelihoods as only median point estimator is used
        /// </summary>
        public ISampleMap<Dictionary<Genotype, double>> GetCopyNumbersLikelihoods(ISampleMap<CanvasSegment> canvasSegments, ISampleMap<SampleMetrics> samplesInfo,
            ISampleMap<ICopyNumberModel> copyNumberModel, int numberOfTrimmedBins)
        {
            var genotypes = Enumerable.Range(0, _maximumCopyNumber).Select(Genotype.Create).ToList();
            const double maxCoverageMultiplier = 3.0;
            var singleSampleLikelihoods = new SampleMap<Dictionary<Genotype, double>>();

            foreach (var sampleId in canvasSegments.SampleIds)
            {
                var density = new Dictionary<Genotype, double>();

                foreach (var genotypeCopyNumber in genotypes)
                {
                    double cvg = Math.Min(canvasSegments[sampleId].TruncatedMedianCount(numberOfTrimmedBins),
                                samplesInfo[sampleId].MeanCoverage * maxCoverageMultiplier);
                    // In case we run into out-of-range trouble again (CANV-694), print details
                    {
                        int intcvg = Convert.ToInt32(cvg);
                        int coverageBound = copyNumberModel[sampleId].GetCoverageBound();
                        double truncatedDepth = canvasSegments[sampleId].TruncatedMedianCount(numberOfTrimmedBins);
                        double meanTimesThree = samplesInfo[sampleId].MeanCoverage * maxCoverageMultiplier;
                        int maxAllowedCN = copyNumberModel[sampleId].GetMaxCopyNumber();
                        if ( intcvg >= coverageBound || genotypeCopyNumber.TotalCopyNumber > maxAllowedCN)
                        {
                            throw new ArgumentException(
                                $"Tried to look up bad depth or CN for {sampleId}: depth {intcvg} CN {genotypeCopyNumber.TotalCopyNumber}" +
                                $" where max handled values are {coverageBound} and {maxAllowedCN} respectively;" +
                                $" original depth was {truncatedDepth}, mean * 3 was {meanTimesThree};" +
                                $" segment {canvasSegments[sampleId].Chr}:{canvasSegments[sampleId].Begin}-{canvasSegments[sampleId].End}");

                        }
                    }
                    double currentLikelihood =
                        copyNumberModel[sampleId].GetTotalCopyNumberLikelihoods(cvg, genotypeCopyNumber);
                    currentLikelihood = Double.IsNaN(currentLikelihood) || Double.IsInfinity(currentLikelihood)
                        ? 0
                        : currentLikelihood;
                    density[genotypeCopyNumber] = currentLikelihood;
                }
                singleSampleLikelihoods.Add(sampleId, density);
            }
            return singleSampleLikelihoods;
        }
    }
}