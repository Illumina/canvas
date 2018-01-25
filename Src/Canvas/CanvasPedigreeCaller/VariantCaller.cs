using System;
using System.Collections.Generic;
using System.Linq;
using CanvasCommon;
using Illumina.Common;
using Isas.Framework.DataTypes;
using Isas.Framework.DataTypes.Maps;
[assembly: System.Runtime.CompilerServices.InternalsVisibleTo("CanvasTest")]

namespace CanvasPedigreeCaller
{
    internal class VariantCaller : IVariantCaller
    {
        private readonly CopyNumberLikelihoodCalculator _copyNumberLikelihoodCalculator;
        private readonly PedigreeCallerParameters _callerParameters;
        private readonly int _qualityFilterThreshold;
        private readonly Dictionary<int, List<PhasedGenotype>> _genotypes;

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
        private static Dictionary<int, List<PhasedGenotype>> GenerateGenotypeCombinations(int numberOfCnStates)
        {
            var genotypes = new Dictionary<int, List<PhasedGenotype>>();
            for (int cn = 0; cn < numberOfCnStates; cn++)
            {
                genotypes[cn] = new List<PhasedGenotype>();
                for (int gt = 0; gt <= cn; gt++)
                {
                    genotypes[cn].Add(new PhasedGenotype(gt, cn - gt));
                }
            }
            return genotypes;
        }

        /// <summary>
        /// Derives metrics from b-allele counts within each segment and determines whereas to use them for calculating MCC
        /// </summary>
        /// <param name="canvasSegments"></param>
        /// <returns></returns>
        private bool UseMafInformation(ISampleMap<CanvasSegment> canvasSegments)
        {
            var alleles = canvasSegments.Values.Select(segments => segments.Balleles?.TotalCoverage);
            var alleleCounts = alleles.Select(allele => allele?.Count ?? 0).ToList();
            bool lowAlleleCounts = alleleCounts.Select(x => x < _callerParameters.DefaultReadCountsThreshold).Any(c => c);
            return lowAlleleCounts;
        }

        private void EstimateQScores(ISampleMap<CanvasSegment> canvasSegments, ISampleMap<SampleMetrics> pedigreeMembersInfo,
            PedigreeInfo pedigreeInfo, ISampleMap<Dictionary<Genotype, double>> SingleSampleLikelihoods, JointLikelihoods copyNumberLikelihoods, ISampleMap<Genotype> copyNumbers)
        {
            foreach (var sampleId in canvasSegments.SampleIds)
            {
                canvasSegments[sampleId].QScore = GetSingleSampleQualityScore(SingleSampleLikelihoods[sampleId], copyNumbers[sampleId]);
                canvasSegments[sampleId].CopyNumber = copyNumbers[sampleId].TotalCopyNumber;
                if (canvasSegments[sampleId].QScore < _qualityFilterThreshold)
                    canvasSegments[sampleId].Filter = CanvasFilter.Create(new[] {$"q{_qualityFilterThreshold}"});
            }
            if (pedigreeInfo.HasFullPedigree())
                SetDenovoQualityScores(canvasSegments, pedigreeMembersInfo, pedigreeInfo.ParentsIds, pedigreeInfo.OffspringIds, copyNumberLikelihoods, copyNumbers);
        }

        private double GetSingleSampleQualityScore(Dictionary<Genotype, double> copyNumbersLikelihoods, Genotype cnState)
        {
            double normalizationConstant = copyNumbersLikelihoods.Select(ll => ll.Value).Sum();
            double qscore = -10.0 * Math.Log10((normalizationConstant - copyNumbersLikelihoods[cnState]) / normalizationConstant);
            if (Double.IsInfinity(qscore) | qscore > _callerParameters.MaxQscore)
                qscore = _callerParameters.MaxQscore;
            return qscore;
        }

        private void SetDenovoQualityScores(ISampleMap<CanvasSegment> canvasSegments, ISampleMap<SampleMetrics> samplesInfo, List<SampleId> parentIDs, List<SampleId> offspringIDs,
            JointLikelihoods copyNumbersLikelihoods, ISampleMap<Genotype> copyNumbers)
        {
            var pedigreeCopyNumbers = copyNumbers.WhereSampleIds(sampleId =>
                parentIDs.Contains(sampleId) || offspringIDs.Contains(sampleId));
            foreach (var probandId in offspringIDs)
            {
                // targeted proband is REF
                if (IsReferenceVariant(canvasSegments, samplesInfo, probandId))
                    continue;
                // common variant
                if (IsCommonCnv(canvasSegments, samplesInfo, parentIDs, probandId, _callerParameters.MaximumCopyNumber))
                    continue;
                // other offsprings are ALT
                if (!offspringIDs.Except(probandId.ToEnumerable()).All(id => IsReferenceVariant(canvasSegments, samplesInfo, id)))
                    continue;
                // not all q-scores are above the threshold
                if (parentIDs.Concat(probandId).Any(id => !IsPassVariant(canvasSegments, id)))
                    continue;

                double deNovoQualityScore = GetConditionalDeNovoQualityScore(copyNumbersLikelihoods, probandId, pedigreeCopyNumbers);
                if (Double.IsInfinity(deNovoQualityScore) | deNovoQualityScore > _callerParameters.MaxQscore)
                    deNovoQualityScore = _callerParameters.MaxQscore;
                canvasSegments[probandId].DqScore = deNovoQualityScore;
            }
        }

        private double GetConditionalDeNovoQualityScore(JointLikelihoods jointLikelihoods, SampleId probandId, ISampleMap<Genotype> copyNumberGenotypes)
        {
            const double q60 = 0.000001;
            double marginalAltLikelihood = jointLikelihoods.GetMarginalLikelihood(new KeyValuePair<SampleId, Genotype>(probandId, copyNumberGenotypes[probandId]));
            double marginalRefLikelihood = jointLikelihoods.GetMarginalLikelihood(new KeyValuePair<SampleId, Genotype>(probandId, Genotype.Create(2)));
            double normalization = marginalAltLikelihood + marginalRefLikelihood;
            double denovoProbability = jointLikelihoods.GetJointLikelihood(copyNumberGenotypes) / marginalAltLikelihood * marginalRefLikelihood / normalization;
            return -10.0 * Math.Log10(Math.Max(denovoProbability, q60));
        }


        private bool IsPassVariant(ISampleMap<CanvasSegment> canvasSegments, SampleId sampleId)
        {
            return canvasSegments[sampleId].QScore >= _qualityFilterThreshold;
        }

        public static bool IsCommonCnv(ISampleMap<CanvasSegment> canvasSegments, ISampleMap<SampleMetrics> samplesInfo, List<SampleId> parentIDs, SampleId probandId, int maximumCopyNumber)
        {
            int parent1CopyNumber = GetCnState(canvasSegments, parentIDs.First(), maximumCopyNumber);
            int parent2CopyNumber = GetCnState(canvasSegments, parentIDs.Last(), maximumCopyNumber);
            int probandCopyNumber = GetCnState(canvasSegments, probandId, maximumCopyNumber);
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
            ISampleMap<ICopyNumberModel> copyNumberModel, PedigreeInfo pedigreeInfo)
        {
            var singleSampleLikelihoods = _copyNumberLikelihoodCalculator.GetCopyNumbersLikelihoods(canvasSegments, samplesInfo, copyNumberModel);
            
            (var pedigreeCopyNumbers, var pedigreeLikelihoods) = GetPedigreeCopyNumbers(pedigreeInfo, singleSampleLikelihoods);

            var nonPedigreeCopyNumbers = CanvasPedigreeCaller.GetNonPedigreeCopyNumbers(canvasSegments, pedigreeInfo, singleSampleLikelihoods);

            var mergedCopyNumbers = pedigreeCopyNumbers.Concat(nonPedigreeCopyNumbers).OrderBy(canvasSegments.SampleIds);

            EstimateQScores(canvasSegments, samplesInfo, pedigreeInfo, singleSampleLikelihoods, pedigreeLikelihoods, mergedCopyNumbers);

            // TODO: this will be integrated with GetCopyNumbers* on a model level as a part of https://jira.illumina.com/browse/CANV-404
            if (!UseMafInformation(canvasSegments) && pedigreeInfo.HasFullPedigree())
                AssignMccWithPedigreeInfo(canvasSegments, copyNumberModel, pedigreeInfo);
            if (!UseMafInformation(canvasSegments) && pedigreeInfo.HasOther())
                AssignMccNoPedigreeInfo(canvasSegments.Where(segment=> pedigreeInfo.OtherIds.Contains(segment.SampleId)).ToSampleMap(), copyNumberModel, _genotypes);
        }
        
        /// <summary>
        /// Calculates maximal likelihood for segments with SNV allele counts given CopyNumber. Updated MajorChromosomeCount.
        /// </summary>   
        private void AssignMccNoPedigreeInfo(ISampleMap<CanvasSegment> canvasSegments,
            ISampleMap<ICopyNumberModel> model, Dictionary<int, List<PhasedGenotype>> genotypes)
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
                double gqscore = GetGtLikelihoodScore(canvasSegments[sampleId].Balleles, genotypeset, ref selectedGtState, model[sampleId]);
                canvasSegments[sampleId].MajorChromosomeCountScore = gqscore;
                if (selectedGtState.HasValue)
                    canvasSegments[sampleId].MajorChromosomeCount =
                        Math.Max(genotypeset[selectedGtState.Value].CopyNumberA,
                            genotypeset[selectedGtState.Value].CopyNumberB);
            }
        }

        /// <summary>
        /// Calculates maximal likelihood for genotypes given a copy number call. Updated MajorChromosomeCount.
        /// </summary>
        private void AssignMccWithPedigreeInfo(ISampleMap<CanvasSegment> canvasSegments, ISampleMap<ICopyNumberModel> model, PedigreeInfo pedigreeInfo)
        {
            double maximalLikelihood = Double.MinValue;
            int parent1CopyNumber = canvasSegments[pedigreeInfo.ParentsIds.First()].CopyNumber;
            int parent2CopyNumber = canvasSegments[pedigreeInfo.ParentsIds.Last()].CopyNumber;

            foreach (var parent1GtStates in _genotypes[parent1CopyNumber])
            {
                foreach (var parent2GtStates in _genotypes[parent2CopyNumber])
                {
                    var bestChildGtStates = new List<PhasedGenotype>();
                    double currentLikelihood = 1;
                    foreach (SampleId child in pedigreeInfo.OffspringIds)
                    {
                        int childCopyNumber = canvasSegments[child].CopyNumber;
                        bool isInheritedCnv = !canvasSegments[child].DqScore.HasValue;
                        double bestLikelihood = Double.MinValue;
                        PhasedGenotype bestGtState = null;
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
                        AssignMcc(canvasSegments[pedigreeInfo.ParentsIds.First()], model[pedigreeInfo.ParentsIds.First()], parent1GtStates, parent1CopyNumber);
                        AssignMcc(canvasSegments[pedigreeInfo.ParentsIds.Last()], model[pedigreeInfo.ParentsIds.Last()], parent2GtStates, parent2CopyNumber);
                        var counter = 0;
                        foreach (SampleId child in pedigreeInfo.OffspringIds)
                        {
                            if (bestChildGtStates[counter] == null) continue;
                            int childCopyNumber = canvasSegments[child].CopyNumber;
                            AssignMcc(canvasSegments[child], model[child], bestChildGtStates[counter], childCopyNumber);
                            counter++;
                        }
                    }
                }
            }
        }

        private double GetProbandLikelihood(ICopyNumberModel copyNumberModel, int childCopyNumber, PhasedGenotype parent1GtStates, PhasedGenotype parent2GtStates, bool isInheritedCnv, CanvasSegment canvasSegment,
            double bestLikelihood, ref PhasedGenotype bestGtState)
        {
            foreach (var childGtState in _genotypes[childCopyNumber])
            {
                double currentChildLikelihood;
                if (IsGtPedigreeConsistent(parent1GtStates, childGtState) &&
                    IsGtPedigreeConsistent(parent2GtStates, childGtState)
                    && isInheritedCnv)
                    currentChildLikelihood = copyNumberModel.GetGenotypeLikelihood(canvasSegment.Balleles, childGtState);
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

        private bool IsGtPedigreeConsistent(PhasedGenotype parentGtStates, PhasedGenotype childGtStates)
        {
            if (parentGtStates.CopyNumberA == childGtStates.CopyNumberA || parentGtStates.CopyNumberB == childGtStates.CopyNumberA ||
                parentGtStates.CopyNumberA == childGtStates.CopyNumberB || parentGtStates.CopyNumberB == childGtStates.CopyNumberB)
                return true;
            return false;
        }
        
        private void AssignMcc(CanvasSegment canvasSegment, ICopyNumberModel copyNumberModel,
            PhasedGenotype gtStates, int copyNumber)
        {
            const int diploidCopyNumber = 2;
            const int haploidCopyNumber = 1;
            if (copyNumber > diploidCopyNumber)
            {

                canvasSegment.MajorChromosomeCount =
                    Math.Max(gtStates.CopyNumberA, gtStates.CopyNumberB);
                int? selectedGtState = _genotypes[copyNumber].IndexOf(gtStates);
                canvasSegment.MajorChromosomeCountScore = GetGtLikelihoodScore(canvasSegment.Balleles, _genotypes[copyNumber], ref selectedGtState, copyNumberModel);
                    GetGtLikelihoodScore(canvasSegment.Balleles, _genotypes[copyNumber], ref selectedGtState, copyNumberModel);
            }
            else
            {
                canvasSegment.MajorChromosomeCount = copyNumber == diploidCopyNumber
                    ? haploidCopyNumber : copyNumber;
                canvasSegment.MajorChromosomeCountScore = null;
            }
        }

        private double GetGtLikelihoodScore(Balleles gtObservedCounts, List<PhasedGenotype> gtModelCounts, ref int? selectedGtState, ICopyNumberModel copyNumberModel)
        {
            const int maxGQscore = 60;
            var gtLikelihoods = Enumerable.Repeat(0.0, gtModelCounts.Count).ToList();
            var gtModelCounter = 0;
            foreach (var gtModelCount in gtModelCounts)
            {
                gtLikelihoods[gtModelCounter] = copyNumberModel.GetGenotypeLikelihood(gtObservedCounts, gtModelCount);
                gtModelCounter++;
            }
            if (!selectedGtState.HasValue)
                selectedGtState = gtLikelihoods.IndexOf(gtLikelihoods.Max());
            double normalizationConstant = gtLikelihoods.Sum();
            double gqscore = -10.0 * Math.Log10((normalizationConstant - gtLikelihoods[selectedGtState.Value]) / normalizationConstant);
            if (Double.IsInfinity(gqscore) | gqscore > maxGQscore)
                gqscore = maxGQscore;
            return Double.IsNaN(gqscore) || Double.IsInfinity(gqscore) ? 0 : gqscore;
        }

        private static double GetCurrentGtLikelihood(ICopyNumberModel copyNumberModel, CanvasSegment canvasSegment, PhasedGenotype gtStates)
        {
            return copyNumberModel.GetGenotypeLikelihood(canvasSegment.Balleles, gtStates);
        }

        private (ISampleMap<Genotype> copyNumbersGenotypes, JointLikelihoods jointLikelihood) GetPedigreeCopyNumbers(PedigreeInfo pedigreeInfo, ISampleMap<Dictionary<Genotype, double>> copyNumbersLikelihoods)
        {
            int nHighestLikelihoodGenotypes = pedigreeInfo != null && pedigreeInfo.OffspringIds.Count >= 2 ? 3 : _callerParameters.MaximumCopyNumber;
            copyNumbersLikelihoods = copyNumbersLikelihoods.SelectValues(l => l.OrderByDescending(kvp => kvp.Value).Take(nHighestLikelihoodGenotypes).ToDictionary());

            var sampleCopyNumbersGenotypes = new SampleMap<Genotype>();
            var jointLikelihood = new JointLikelihoods();
            if (!pedigreeInfo.HasFullPedigree())
                return (sampleCopyNumbersGenotypes, jointLikelihood);
            foreach (var copyNumberParent1 in copyNumbersLikelihoods[pedigreeInfo.ParentsIds.First()])
            {
                foreach (var copyNumberParent2 in copyNumbersLikelihoods[pedigreeInfo.ParentsIds.Last()])
                {
                    foreach (var offspringGtStates in pedigreeInfo.OffspringPhasedGenotypes)
                    {
                        if (!pedigreeInfo.OffspringIds.All(id => copyNumbersLikelihoods[id].ContainsKey(
                            Genotype.Create(Math.Min(offspringGtStates[pedigreeInfo.OffspringIds.IndexOf(id)].PhasedGenotype.CopyNumberA + offspringGtStates[pedigreeInfo.OffspringIds.IndexOf(id)].PhasedGenotype.CopyNumberB,
                                _callerParameters.MaximumCopyNumber - 1)))))
                        {
                            continue;
                        }
                        double currentLikelihood = copyNumberParent1.Value * copyNumberParent2.Value;
                        var totalCopyNumberGenotypes = new List<Genotype>();
                        for (var counter = 0; counter < pedigreeInfo.OffspringIds.Count; counter++)
                        {
                            var child = pedigreeInfo.OffspringIds[counter];
                            var copyNumberGenotypeChild = Genotype.Create(Math.Min(offspringGtStates[counter].PhasedGenotype.CopyNumberA + offspringGtStates[counter].PhasedGenotype.CopyNumberB,
                                _callerParameters.MaximumCopyNumber - 1));
                            totalCopyNumberGenotypes.Add(copyNumberGenotypeChild);
                            currentLikelihood *= pedigreeInfo.TransitionMatrix[copyNumberParent1.Key.TotalCopyNumber][offspringGtStates[counter].PhasedGenotype.CopyNumberA] *
                                                 pedigreeInfo.TransitionMatrix[copyNumberParent2.Key.TotalCopyNumber][offspringGtStates[counter].PhasedGenotype.CopyNumberB] *
                                                 copyNumbersLikelihoods[child][copyNumberGenotypeChild];
                        }

                        currentLikelihood = Double.IsNaN(currentLikelihood) || Double.IsInfinity(currentLikelihood) ? 0 : currentLikelihood;
                        var genotypesInPedigree = new SampleMap<Genotype>
                        {
                            {pedigreeInfo.ParentsIds.First(), copyNumberParent1.Key},
                            {pedigreeInfo.ParentsIds.Last(), copyNumberParent2.Key}
                        };
                        pedigreeInfo.OffspringIds.Zip(totalCopyNumberGenotypes).ForEach(sampleIdGenotypeKvp => genotypesInPedigree.Add(sampleIdGenotypeKvp.Item1, sampleIdGenotypeKvp.Item2));
                        genotypesInPedigree = genotypesInPedigree.OrderBy(pedigreeInfo.AllSampleIds);
                        jointLikelihood.AddJointLikelihood(genotypesInPedigree, currentLikelihood);

                        if (currentLikelihood > jointLikelihood.MaximalLikelihood)
                        {
                            jointLikelihood.MaximalLikelihood = currentLikelihood;
                            sampleCopyNumbersGenotypes = genotypesInPedigree;
                        }
                    }
                }
            }
            if (sampleCopyNumbersGenotypes.Empty())
                throw new IlluminaException("Maximal likelihood was not found");
            return (sampleCopyNumbersGenotypes, jointLikelihood);
        }
    }
}