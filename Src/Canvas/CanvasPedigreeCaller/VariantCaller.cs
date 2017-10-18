using System;
using System.Collections.Generic;
using System.Linq;
using CanvasCommon;
using Illumina.Common;
using Isas.Framework.DataTypes;
using Isas.Framework.DataTypes.Maps;

namespace CanvasPedigreeCaller
{
    class VariantCaller
    {
        protected double MedianCoverageThreshold = 4;
        private readonly CopyNumberLikelihoodCalculator _copyNumberLikelihoodCalculator;
        private readonly PedigreeCallerParameters _callerParameters;
        private readonly int _qualityFilterThreshold;
        private readonly Dictionary<int, List<Genotype>> _genotypes;

        public VariantCaller(CopyNumberLikelihoodCalculator copyNumberLikelihoodCalculator, PedigreeCallerParameters callerParameters, int qualityFilterThreshold)
        {
            _copyNumberLikelihoodCalculator = copyNumberLikelihoodCalculator;
            _callerParameters = callerParameters;
            _qualityFilterThreshold = qualityFilterThreshold;
            _genotypes = GenerateGenotypeCombinations(callerParameters.MaximumCopyNumber);
        }
        
        /// <summary>
        /// Generate all possible copy number genotype combinations with the maximal number of alleles per segment set to maxAlleleNumber.
        /// </summary>
        /// <param name="numberOfCnStates"></param>
        /// <returns> </returns>
        private static Dictionary<int, List<Genotype>> GenerateGenotypeCombinations(int numberOfCnStates)
        {
            var genotypes = new Dictionary<int, List<Genotype>>();
            for (int cn = 0; cn < numberOfCnStates; cn++)
            {
                genotypes[cn] = new List<Genotype>();
                for (int gt = 0; gt <= cn; gt++)
                {
                    genotypes[cn].Add(new Genotype(gt, cn - gt));
                }
            }
            return genotypes;
        }

        /// <summary>
        /// Derives metrics from b-allele counts within each segment and determines whereas to use them for calculating MCC
        /// </summary>
        /// <param name="canvasSegments"></param>
        /// <param name="segmentIndex"></param>
        /// <returns></returns>
        private bool UseMafInformation(ISampleMap<CanvasSegment> canvasSegments)
        {
            var alleles = canvasSegments.Values.Select(segments => segments.Balleles?.TotalCoverage);
            var alleleCounts = alleles.Select(allele => allele?.Count ?? 0).ToList();
            bool lowAlleleCounts = alleleCounts.Select(x => x < _callerParameters.DefaultReadCountsThreshold).Any(c => c == true);
            var coverageCounts = canvasSegments.Values.Select(segments => segments.MedianCount).ToList();
            var isSkewedHetHomRatio = false;
            double alleleDensity = canvasSegments.Values.First().Length /
                                   Math.Max(alleleCounts.Average(), 1.0);
            bool useCnLikelihood = lowAlleleCounts ||
                                   alleleDensity < _callerParameters.DefaultAlleleDensityThreshold ||
                                   alleleCounts.Any(x => x > _callerParameters.DefaultPerSegmentAlleleMaxCounts) ||
                                   coverageCounts.Any(coverage => coverage < MedianCoverageThreshold) ||
                                   isSkewedHetHomRatio;
            // for now only use lowAlleleCounts metric
            return lowAlleleCounts;
        }

        private void EstimateQScores(ISampleMap<CanvasSegment> canvasSegments, ISampleMap<SampleMetrics> pedigreeMembersInfo,
            PedigreeInfo pedigreeInfo, CopyNumbersLikelihoods copyNumberLikelihoods, ISampleMap<int> copyNumbers)
        {
            foreach (var sampleId in canvasSegments.SampleIds)
            {
                canvasSegments[sampleId].QScore = GetSingleSampleQualityScore(copyNumberLikelihoods.SingleSampleLikelihoods[sampleId], copyNumbers[sampleId]);
                canvasSegments[sampleId].CopyNumber = copyNumbers[sampleId];
                if (canvasSegments[sampleId].QScore < _qualityFilterThreshold)
                    canvasSegments[sampleId].Filter = $"q{_qualityFilterThreshold}";
            }
            if (pedigreeInfo != null)
                SetDenovoQualityScores(canvasSegments, pedigreeMembersInfo, pedigreeInfo.ParentsIds, pedigreeInfo.OffspringIds, copyNumberLikelihoods);
        }

        private double GetSingleSampleQualityScore(Dictionary<int, double> copyNumbersLikelihoods, int cnState)
        {
            double normalizationConstant = copyNumbersLikelihoods.Select(ll => ll.Value).Sum();
            double qscore = -10.0 * Math.Log10((normalizationConstant - copyNumbersLikelihoods[cnState]) / normalizationConstant);
            if (Double.IsInfinity(qscore) | qscore > _callerParameters.MaxQscore)
                qscore = _callerParameters.MaxQscore;
            return qscore;
        }

        private void SetDenovoQualityScores(ISampleMap<CanvasSegment> canvasSegments, ISampleMap<SampleMetrics> samplesInfo, List<SampleId> parentIDs, List<SampleId> offspringIDs,
            CopyNumbersLikelihoods copyNumberLikelihoods)
        {
            foreach (var probandId in offspringIDs)
            {
                // targeted proband is REF
                if (IsReferenceVariant(canvasSegments, samplesInfo, probandId))
                    continue;
                // common variant
                if (IsCommonCnv(canvasSegments, samplesInfo, parentIDs, probandId))
                    continue;
                // other offsprings are ALT
                if (!offspringIDs.Except(probandId.ToEnumerable()).All(id => IsReferenceVariant(canvasSegments, samplesInfo, id)))
                    continue;
                // not all q-scores are above the threshold
                if (parentIDs.Concat(probandId).Any(id => !IsPassVariant(canvasSegments, id)))
                    continue;

                double deNovoQualityScore = GetConditionalDeNovoQualityScore(copyNumberLikelihoods, probandId, canvasSegments, parentIDs);
                if (Double.IsInfinity(deNovoQualityScore) | deNovoQualityScore > _callerParameters.MaxQscore)
                    deNovoQualityScore = _callerParameters.MaxQscore;
                canvasSegments[probandId].DqScore = deNovoQualityScore;
            }
        }

        private double GetConditionalDeNovoQualityScore(CopyNumbersLikelihoods copyNumbersLikelihoods, SampleId probandId, ISampleMap<CanvasSegment> canvasSegments,
            List<SampleId> parentIDs)
        {
            var numerator = 0.0;
            var denominator = 0.0;
            const int diploidState = 2;
            var probandCopyNumber = canvasSegments[probandId].CopyNumber;
            var probandMarginalProbabilities = copyNumbersLikelihoods.GetMarginalProbability(_callerParameters.MaximumCopyNumber, probandId);
            double normalization = probandMarginalProbabilities[probandCopyNumber] + probandMarginalProbabilities[diploidState];
            double probandMarginalAlt = probandMarginalProbabilities[probandCopyNumber] / normalization;

            foreach (var copyNumberIndex in copyNumbersLikelihoods.CopyNumbers.Where(copyNumbers => copyNumbers[probandId] == canvasSegments[probandId].CopyNumber))
            {
                double holder = copyNumbersLikelihoods.GetJointLikelihood(copyNumberIndex);
                denominator += holder;
                if (copyNumberIndex[parentIDs.First()] == diploidState && copyNumberIndex[parentIDs.Last()] == diploidState)
                    numerator += holder;
            }

            const double q60 = 0.000001;
            double denovoProbability = (1 - numerator / denominator) * (1 - probandMarginalAlt);
            return -10.0 * Math.Log10(Math.Max(denovoProbability, q60));
        }

        private bool IsPassVariant(ISampleMap<CanvasSegment> canvasSegments, SampleId sampleId)
        {
            return canvasSegments[sampleId].QScore > _qualityFilterThreshold;
        }

        private bool IsCommonCnv(ISampleMap<CanvasSegment> canvasSegments, ISampleMap<SampleMetrics> samplesInfo, List<SampleId> parentIDs, SampleId probandId)
        {
            int parent1CopyNumber = GetCnState(canvasSegments, parentIDs.First(), _callerParameters.MaximumCopyNumber);
            int parent2CopyNumber = GetCnState(canvasSegments, parentIDs.Last(), _callerParameters.MaximumCopyNumber);
            int probandCopyNumber = GetCnState(canvasSegments, probandId, _callerParameters.MaximumCopyNumber);
            var parent1Genotypes = CanvasPedigreeCaller.GenerateCnAlleles(parent1CopyNumber);
            var parent2Genotypes = CanvasPedigreeCaller.GenerateCnAlleles(parent2CopyNumber);
            var probandGenotypes = CanvasPedigreeCaller.GenerateCnAlleles(probandCopyNumber);
            var parent1Segment = canvasSegments[parentIDs.First()];
            var parent2Segment = canvasSegments[parentIDs.Last()];
            var probandSegment = canvasSegments[probandId];
            int parent1Ploidy = samplesInfo[probandId].GetPloidy(parent1Segment);
            int parent2Ploidy = samplesInfo[probandId].GetPloidy(parent2Segment);
            int probandPloidy = samplesInfo[probandId].GetPloidy(probandSegment);
            bool isCommoCnv = parent1Genotypes.Intersect(probandGenotypes).Any() && parent1Ploidy == probandPloidy ||
                              parent2Genotypes.Intersect(probandGenotypes).Any() && parent2Ploidy == probandPloidy;
            return isCommoCnv;
        }

        private bool IsReferenceVariant(ISampleMap<CanvasSegment> canvasSegments, ISampleMap<SampleMetrics> samplesInfo, SampleId sampleId)
        {
            var segment = canvasSegments[sampleId];
            return GetCnState(canvasSegments, sampleId, _callerParameters.MaximumCopyNumber) == samplesInfo[sampleId].GetPloidy(segment);
        }

        private static int GetCnState(ISampleMap<CanvasSegment> canvasSegmentsSet, SampleId sampleId, int maximumCopyNumber)
        {
            return Math.Min(canvasSegmentsSet[sampleId].CopyNumber, maximumCopyNumber - 1);
        }

        /// <summary>
        /// Identify variant with the highest likelihood at a given setPosition and assign relevant scores
        /// </summary>
        public void CallVariant(ISampleMap<CanvasSegment> canvasSegments, ISampleMap<SampleMetrics> samplesInfo,
            ISampleMap<CopyNumberModel> copyNumberModel,
            PedigreeInfo pedigreeInfo)
        {
            var copyNumbersLikelihoods = _copyNumberLikelihoodCalculator.GetCopyNumbersLikelihoods(canvasSegments, samplesInfo, copyNumberModel);

            var copyNumbers = pedigreeInfo != null
                ? GetCopyNumbersWithPedigreeInfo(canvasSegments, copyNumberModel, pedigreeInfo, copyNumbersLikelihoods)
                : CanvasPedigreeCaller.GetCopyNumbersNoPedigreeInfo(canvasSegments, copyNumbersLikelihoods);

            EstimateQScores(canvasSegments, samplesInfo, pedigreeInfo, copyNumbersLikelihoods, copyNumbers);

            // TODO: this will be integrated with GetCopyNumbers* on a model level as a part of https://jira.illumina.com/browse/CANV-404
            if (!UseMafInformation(canvasSegments) && pedigreeInfo != null)
                AssignMccWithPedigreeInfo(canvasSegments, copyNumberModel, pedigreeInfo);
            if (!UseMafInformation(canvasSegments) && pedigreeInfo == null)
                AssignMccNoPedigreeInfo(canvasSegments, copyNumberModel, _genotypes);
        }
        
        /// <summary>
        /// Calculates maximal likelihood for segments with SNV allele counts given CopyNumber. Updated MajorChromosomeCount.
        /// </summary>   
        private void AssignMccNoPedigreeInfo(ISampleMap<CanvasSegment> canvasSegments,
            ISampleMap<CopyNumberModel> model, Dictionary<int, List<Genotype>> genotypes)
        {
            foreach (var sampleId in canvasSegments.SampleIds)
            {
                int copyNumber = canvasSegments[sampleId].CopyNumber;
                if (copyNumber > 2)
                {
                    canvasSegments[sampleId].MajorChromosomeCount = copyNumber == 2 ? 1 : copyNumber;
                    return;
                }
                var genotypeset = genotypes[copyNumber];
                int? selectedGtState = null;
                double gqscore = model[sampleId].GetGtLikelihoodScore(canvasSegments[sampleId].Balleles.GetAlleleCounts(),
                    genotypeset, ref selectedGtState);
                canvasSegments[sampleId].MajorChromosomeCountScore = gqscore;
                if (selectedGtState.HasValue)
                    canvasSegments[sampleId].MajorChromosomeCount =
                        Math.Max(genotypeset[selectedGtState.Value].CountsA,
                            genotypeset[selectedGtState.Value].CountsB);
            }
        }

        /// <summary>
        /// Calculates maximal likelihood for genotypes given a copy number call. Updated MajorChromosomeCount.
        /// </summary>
        private void AssignMccWithPedigreeInfo(ISampleMap<CanvasSegment> canvasSegments, ISampleMap<CopyNumberModel> model, PedigreeInfo pedigreeInfo)
        {
            double maximalLikelihood = Double.MinValue;
            int parent1CopyNumber = canvasSegments[pedigreeInfo.ParentsIds.First()].CopyNumber;
            int parent2CopyNumber = canvasSegments[pedigreeInfo.ParentsIds.Last()].CopyNumber;

            foreach (var parent1GtStates in _genotypes[parent1CopyNumber])
            {
                foreach (var parent2GtStates in _genotypes[parent2CopyNumber])
                {
                    var bestChildGtStates = new List<Genotype>();
                    double currentLikelihood = 1;
                    foreach (SampleId child in pedigreeInfo.OffspringIds)
                    {
                        int childCopyNumber = canvasSegments[child].CopyNumber;
                        bool isInheritedCnv = !canvasSegments[child].DqScore.HasValue;
                        double bestLikelihood = Double.MinValue;
                        Genotype bestGtState = null;
                        bestLikelihood = GetProbandLikelihood(model[child], childCopyNumber,
                            parent1GtStates, parent2GtStates, isInheritedCnv, canvasSegments[child], bestLikelihood, ref bestGtState);
                        bestChildGtStates.Add(bestGtState);
                        currentLikelihood *= bestLikelihood;
                    }
                    currentLikelihood *= GetCurrentGtLikelihood(model[pedigreeInfo.ParentsIds.First()], canvasSegments[pedigreeInfo.ParentsIds.First()], parent1GtStates) *
                                         GetCurrentGtLikelihood(model[pedigreeInfo.ParentsIds.Last()], canvasSegments[pedigreeInfo.ParentsIds.Last()], parent2GtStates);

                    currentLikelihood = Double.IsNaN(currentLikelihood) || Double.IsInfinity(currentLikelihood)
                        ? 0
                        : currentLikelihood;

                    if (currentLikelihood > maximalLikelihood)
                    {
                        maximalLikelihood = currentLikelihood;
                        AssignMCC(canvasSegments[pedigreeInfo.ParentsIds.First()], model[pedigreeInfo.ParentsIds.First()], parent1GtStates, parent1CopyNumber);
                        AssignMCC(canvasSegments[pedigreeInfo.ParentsIds.Last()], model[pedigreeInfo.ParentsIds.Last()], parent2GtStates, parent2CopyNumber);
                        var counter = 0;
                        foreach (SampleId child in pedigreeInfo.OffspringIds)
                        {
                            if (bestChildGtStates[counter] == null) continue;
                            int childCopyNumber = canvasSegments[child].CopyNumber;
                            AssignMCC(canvasSegments[child], model[child], bestChildGtStates[counter], childCopyNumber);
                            counter++;
                        }
                    }
                }
            }
        }
        
        private double GetProbandLikelihood(CopyNumberModel copyNumberModel, int childCopyNumber, Genotype parent1GtStates, Genotype parent2GtStates, bool isInheritedCnv, CanvasSegment canvasSegment,
            double bestLikelihood, ref Genotype bestGtState)
        {
            foreach (var childGtState in _genotypes[childCopyNumber])
            {
                double currentChildLikelihood;
                if (IsGtPedigreeConsistent(parent1GtStates, childGtState) &&
                    IsGtPedigreeConsistent(parent2GtStates, childGtState)
                    && isInheritedCnv)
                    currentChildLikelihood = copyNumberModel.GetCurrentGtLikelihood(canvasSegment.Balleles.GetAlleleCounts(), childGtState);
                else
                    continue;
                if (currentChildLikelihood > bestLikelihood)
                {
                    bestLikelihood = currentChildLikelihood;
                    bestGtState = childGtState;
                }
            }
            return bestLikelihood;
        }

        private bool IsGtPedigreeConsistent(Genotype parentGtStates, Genotype childGtStates)
        {
            if (parentGtStates.CountsA == childGtStates.CountsA || parentGtStates.CountsB == childGtStates.CountsA ||
                parentGtStates.CountsA == childGtStates.CountsB || parentGtStates.CountsB == childGtStates.CountsB)
                return true;
            return false;
        }
        
        private void AssignMCC(CanvasSegment canvasSegment, CopyNumberModel copyNumberModel,
            Genotype gtStates, int copyNumber)
        {
            const int diploidCopyNumber = 2;
            const int haploidCopyNumber = 1;
            if (copyNumber > diploidCopyNumber)
            {

                canvasSegment.MajorChromosomeCount =
                    Math.Max(gtStates.CountsA, gtStates.CountsB);
                int? selectedGtState = _genotypes[copyNumber].IndexOf(gtStates);
                canvasSegment.MajorChromosomeCountScore =
                    copyNumberModel.GetGtLikelihoodScore(canvasSegment.Balleles.GetAlleleCounts(), _genotypes[copyNumber], ref selectedGtState);
                copyNumberModel.GetGtLikelihoodScore(canvasSegment.Balleles.GetAlleleCounts(), _genotypes[copyNumber], ref selectedGtState);
            }
            else
            {
                canvasSegment.MajorChromosomeCount = copyNumber == diploidCopyNumber
                    ? haploidCopyNumber : copyNumber;
                canvasSegment.MajorChromosomeCountScore = null;
            }
        }

        private static double GetCurrentGtLikelihood(CopyNumberModel copyNumberModel, CanvasSegment canvasSegment, Genotype gtStates)
        {
            return copyNumberModel.GetCurrentGtLikelihood(canvasSegment.Balleles.GetAlleleCounts(), gtStates);
        }

        /// <summary>
        /// Calculates maximal likelihood for copy numbers. Updated CanvasSegment CopyNumber only. 
        /// </summary>
        private ISampleMap<int> GetCopyNumbersWithPedigreeInfo(ISampleMap<CanvasSegment> segments, ISampleMap<CopyNumberModel> model,
            PedigreeInfo pedigreeInfo, CopyNumbersLikelihoods copyNumbersLikelihoods)
        {
            var sampleCopyNumbers = new SampleMap<int>();
            segments.SampleIds.ForEach(id => sampleCopyNumbers.Add(id, 2));
            var parent1Likelihood = copyNumbersLikelihoods.SingleSampleLikelihoods[pedigreeInfo.ParentsIds.First()];
            var parent2Likelihood = copyNumbersLikelihoods.SingleSampleLikelihoods[pedigreeInfo.ParentsIds.Last()];
            var copyNumbersRange = Enumerable.Range(0, _callerParameters.MaximumCopyNumber).ToList();
            foreach (int copyNumberParent1 in copyNumbersRange)
            {
                foreach (int copyNumberParent2 in copyNumbersRange)
                {
                    foreach (var offspringGtStates in pedigreeInfo.OffspringGenotypes)
                    {
                        double currentLikelihood = parent1Likelihood[copyNumberParent1] * parent2Likelihood[copyNumberParent2];
                        for (var counter = 0; counter < pedigreeInfo.OffspringIds.Count; counter++)
                        {
                            var child = pedigreeInfo.OffspringIds[counter];
                            int copyNumberChild = Math.Min(offspringGtStates[counter].CountsA + offspringGtStates[counter].CountsB,
                                _callerParameters.MaximumCopyNumber - 1);
                            currentLikelihood *= pedigreeInfo.TransitionMatrix[copyNumberParent1][offspringGtStates[counter].CountsA] *
                                                 pedigreeInfo.TransitionMatrix[copyNumberParent2][offspringGtStates[counter].CountsB] *
                                                 copyNumbersLikelihoods.SingleSampleLikelihoods[child][copyNumberChild];
                        }

                        currentLikelihood = Double.IsNaN(currentLikelihood) || Double.IsInfinity(currentLikelihood) ? 0 : currentLikelihood;
                        copyNumbersLikelihoods.SetJointLikelihood(currentLikelihood, GetSampleCopyNumbers(pedigreeInfo, copyNumberParent1, copyNumberParent2, offspringGtStates));

                        if (currentLikelihood > copyNumbersLikelihoods.MaximalLikelihood)
                        {
                            copyNumbersLikelihoods.MaximalLikelihood = currentLikelihood;
                            sampleCopyNumbers = GetSampleCopyNumbers(pedigreeInfo, copyNumberParent1, copyNumberParent2, offspringGtStates);
                        }
                    }
                }
            }
            return sampleCopyNumbers;
        }
        
        private static SampleMap<int> GetSampleCopyNumbers(PedigreeInfo pedigreeInfo, int copyNumberParent1, int copyNumberParent2, List<Genotype> offspringGtStates)
        {
            var sampleCopyNumbers = new SampleMap<int>
            {
                {pedigreeInfo.ParentsIds.First(), copyNumberParent1},
                {pedigreeInfo.ParentsIds.Last(), copyNumberParent2}
            };
            for (int counter = 0; counter < pedigreeInfo.OffspringIds.Count; counter++)
            {
                sampleCopyNumbers.Add(pedigreeInfo.OffspringIds[counter],
                    offspringGtStates[counter].CountsA + offspringGtStates[counter].CountsB);
            }
            return sampleCopyNumbers;
        }
    }
}