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
        public CopyNumbersLikelihoods GetCopyNumbersLikelihoods2(ISampleMap<CanvasSegment> canvasSegments, ISampleMap<SampleMetrics> samplesInfo,
            ISampleMap<ICopyNumberModel> copyNumberModel)
        {
            const double maxCoverageMultiplier = 3.0;
            var singleSampleLikelihoods = new SampleMap<Dictionary<int, double>>();

            foreach (var sampleId in canvasSegments.SampleIds)
            {
                var density = new Dictionary<int, double>();

                foreach (int copyNumber in Enumerable.Range(0, _maximumCopyNumber))
                {
                    double currentLikelihood =
                        copyNumberModel[sampleId].GetTotalCopyNumberLikelihoods(
                            Math.Min(canvasSegments[sampleId].MedianCount,
                                samplesInfo[sampleId].MeanCoverage * maxCoverageMultiplier), Genotype.Create(copyNumber));
                    currentLikelihood = Double.IsNaN(currentLikelihood) || Double.IsInfinity(currentLikelihood)
                        ? 0
                        : currentLikelihood;
                    density[copyNumber] = currentLikelihood;
                }
                singleSampleLikelihoods.Add(sampleId, density);
            }
            return new CopyNumbersLikelihoods(singleSampleLikelihoods, _maximumCopyNumber);
        }
    
    /// <summary>
    /// Calculates maximal likelihood for segments without SNV allele ratios. Updated CanvasSegment CopyNumber only. 
    /// </summary>
    public ISampleMap<Dictionary<Genotype, double>> GetCopyNumbersLikelihoods(ISampleMap<CanvasSegment> canvasSegments, ISampleMap<SampleMetrics> samplesInfo,
            ISampleMap<ICopyNumberModel> copyNumberModel)
        {
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
                            Math.Min(canvasSegments[sampleId].MedianCount,
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