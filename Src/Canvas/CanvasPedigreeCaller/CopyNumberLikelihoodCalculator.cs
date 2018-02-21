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
            ISampleMap<ICopyNumberModel> copyNumberModel)
        {
            const double proportionBins2remove = 0.1;
            var genotypes = Enumerable.Range(0, _maximumCopyNumber).Select(Genotype.Create).ToList();
            const double maxCoverageMultiplier = 3.0;
            var singleSampleLikelihoods = new SampleMap<Dictionary<Genotype, double>>();

            foreach (var sampleId in canvasSegments.SampleIds)
            {
                var density = new Dictionary<Genotype, double>();

                foreach (var genotypeCopyNumber in genotypes)
                {
                    double currentLikelihood =
                        copyNumberModel[sampleId].GetTotalCopyNumberLikelihoods(
                            Math.Min(canvasSegments[sampleId].TruncatedMedianCount(proportionBins2remove),
                                samplesInfo[sampleId].MeanCoverage * maxCoverageMultiplier), genotypeCopyNumber);
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