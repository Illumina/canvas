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

        private void EstimateQScores(ISampleMap<CanvasSegment> canvasSegments, ISampleMap<SampleMetrics> pedigreeMembersInfo,
            PedigreeInfo pedigreeInfo, ISampleMap<Dictionary<Genotype, double>> singleSampleLikelihoods, JointLikelihoods copyNumberLikelihoods, ISampleMap<Genotype> copyNumbers)
        {
            foreach (var sampleId in canvasSegments.SampleIds)
            {
                canvasSegments[sampleId].QScore = GetSingleSampleQualityScore(singleSampleLikelihoods[sampleId], copyNumbers[sampleId]);
                canvasSegments[sampleId].CopyNumber = copyNumbers[sampleId].TotalCopyNumber;
                if (canvasSegments[sampleId].QScore < _qualityFilterThreshold)
                    canvasSegments[sampleId].Filter = CanvasFilter.Create(new[] { $"q{_qualityFilterThreshold}" });
            }
            if (pedigreeInfo.HasFullPedigree())
                SetDenovoQualityScores(canvasSegments, pedigreeMembersInfo, pedigreeInfo.ParentsIds, pedigreeInfo.OffspringIds, copyNumberLikelihoods);
        }

        private double GetSingleSampleQualityScore(Dictionary<Genotype, double> copyNumbersLikelihoods, Genotype cnState)
        {
            double normalizationConstant = copyNumbersLikelihoods.Select(ll => ll.Value).Sum();
            double qscore = -10.0 * Math.Log10((normalizationConstant - copyNumbersLikelihoods[cnState]) / normalizationConstant);
            if (Double.IsInfinity(qscore) | qscore > _callerParameters.MaxQscore)
                qscore = _callerParameters.MaxQscore;
            return qscore;
        }

        /// <summary>
        /// Perform de-novo CNV calling in two steps:
        /// 1. Filter REF variants and common CNVs, this step relies only on total CN calls with associated shortcomings 
        /// 2. Assign de-novo quality based on joint likelihood across pedigree using marginalisation operations  
        /// </summary>
        /// <param name="canvasSegments"></param>
        /// <param name="samplesInfo"></param>
        /// <param name="parentIDs"></param>
        /// <param name="offspringIDs"></param>
        /// <param name="copyNumbersLikelihoods"></param>
        private void SetDenovoQualityScores(ISampleMap<CanvasSegment> canvasSegments, ISampleMap<SampleMetrics> samplesInfo, List<SampleId> parentIDs, List<SampleId> offspringIDs,
            JointLikelihoods copyNumbersLikelihoods)
        {
            foreach (var probandId in offspringIDs)
            {
                // targeted proband is REF
                if (IsReferenceVariant(canvasSegments, samplesInfo, probandId))
                    continue;
                // common variant
                if (IsSharedCnv(canvasSegments, samplesInfo, parentIDs, probandId, _callerParameters.MaximumCopyNumber))
                    continue;
                // other offsprings are ALT
                if (!offspringIDs.Except(probandId.ToEnumerable()).All(id => IsReferenceVariant(canvasSegments, samplesInfo, id)))
                    continue;
                // not all q-scores are above the threshold
                if (parentIDs.Concat(probandId).Any(id => !IsPassVariant(canvasSegments, id)))
                    continue;
                double deNovoQualityScore = GetConditionalDeNovoQualityScore(canvasSegments, copyNumbersLikelihoods, samplesInfo, parentIDs, probandId);

                // adjustment so that denovo quality score threshold is 20 (rather than 10) to match Manta 
                deNovoQualityScore *= 2;

                if (Double.IsInfinity(deNovoQualityScore) | deNovoQualityScore > _callerParameters.MaxQscore)
                    deNovoQualityScore = _callerParameters.MaxQscore;
                canvasSegments[probandId].DqScore = deNovoQualityScore;
            }
        }

        /// <summary>
        /// Assess likelihood of a de-novo variant for copyNumberGenotypes configuration with a Mendelian conflict 
        /// </summary>
        /// <param name="canvasSegments"></param>
        /// <param name="jointLikelihoods"></param>
        /// <param name="parentIDs"></param>
        /// <param name="probandId"></param>
        /// <param name="samplesInfo"></param>
        /// <returns></returns>
        private double GetConditionalDeNovoQualityScore(ISampleMap<CanvasSegment> canvasSegments, JointLikelihoods jointLikelihoods, ISampleMap<SampleMetrics> samplesInfo, List<SampleId> parentIDs, SampleId probandId)
        {
            const double q60 = 0.000001;
            var parent1Ploidy = Genotype.Create(samplesInfo[parentIDs.First()].GetPloidy(canvasSegments[parentIDs.First()]));
            var parent2Ploidy = Genotype.Create(samplesInfo[parentIDs.Last()].GetPloidy(canvasSegments[parentIDs.Last()]));
            int probandPloidy = samplesInfo[probandId].GetPloidy(canvasSegments[probandId]);

            double deNovoGainMarginalLikelihood = jointLikelihoods.GetMarginalGainDeNovoLikelihood(new KeyValuePair<SampleId, Genotype>(probandId, Genotype.Create(probandPloidy)),
                    new KeyValuePair<SampleId, Genotype>(parentIDs.First(), parent1Ploidy), new KeyValuePair<SampleId, Genotype>(parentIDs.Last(), parent2Ploidy));
            double deNovoLossMarginalLikelihood = jointLikelihoods.GetMarginalLossDeNovoLikelihood(new KeyValuePair<SampleId, Genotype>(probandId, Genotype.Create(probandPloidy)),
                    new KeyValuePair<SampleId, Genotype>(parentIDs.First(), parent1Ploidy), new KeyValuePair<SampleId, Genotype>(parentIDs.Last(), parent2Ploidy));
            double denovoProbability = canvasSegments[probandId].CopyNumber > probandPloidy ?
                1 - deNovoGainMarginalLikelihood / (jointLikelihoods.TotalMarginalLikelihood - deNovoLossMarginalLikelihood) :
                1 - deNovoLossMarginalLikelihood / (jointLikelihoods.TotalMarginalLikelihood - deNovoGainMarginalLikelihood);
            // likelihood of proband genotype != ALT given "copyNumberGenotypes" configuration in pedigree with Mendelian conflict 
            return -10.0 * Math.Log10(Math.Max(denovoProbability, q60));
        }


        private bool IsPassVariant(ISampleMap<CanvasSegment> canvasSegments, SampleId sampleId)
        {
            return canvasSegments[sampleId].QScore >= _qualityFilterThreshold;
        }

        /// <summary>
        /// identify common variants using total CN calls within a pedigree obtained with coverage information only 
        /// </summary>
        /// <param name="canvasSegments"></param>
        /// <param name="samplesInfo"></param>
        /// <param name="parentIDs"></param>
        /// <param name="probandId"></param>
        /// <param name="maximumCopyNumber"></param>
        /// <returns></returns>
        public static bool IsSharedCnv(ISampleMap<CanvasSegment> canvasSegments, ISampleMap<SampleMetrics> samplesInfo, List<SampleId> parentIDs, SampleId probandId, int maximumCopyNumber)
        {
            int parent1CopyNumber = GetCnState(canvasSegments, parentIDs.First(), maximumCopyNumber);
            int parent2CopyNumber = GetCnState(canvasSegments, parentIDs.Last(), maximumCopyNumber);
            int probandCopyNumber = GetCnState(canvasSegments, probandId, maximumCopyNumber);
            var parent1Segment = canvasSegments[parentIDs.First()];
            var parent2Segment = canvasSegments[parentIDs.Last()];
            var probandSegment = canvasSegments[probandId];
            int parent1Ploidy = samplesInfo[parentIDs.First()].GetPloidy(parent1Segment);
            int parent2Ploidy = samplesInfo[parentIDs.Last()].GetPloidy(parent2Segment);
            int probandPloidy = samplesInfo[probandId].GetPloidy(probandSegment);
            // Use the following logic: if the proband has fewer copies than expected (from ploidy) but both parents have at least the expected number of copies OR the 
            // proband has more copies than expected but both parents have no more than the expected number of copies, 
            // then it is not a 'common CNV' (i.e.it could be de novo); otherwise, it is common
            return !(parent1CopyNumber <= parent1Ploidy && parent2CopyNumber <= parent2Ploidy && probandCopyNumber > probandPloidy ||
                parent1CopyNumber >= parent1Ploidy && parent2CopyNumber >= parent2Ploidy && probandCopyNumber < probandPloidy);
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
            if (CanvasPedigreeCaller.UseAlleleCountsInformation(canvasSegments, _callerParameters.MinAlleleCountsThreshold, _callerParameters.MinAlleleNumberInSegment) &&
                pedigreeInfo.HasFullPedigree())
                AssignMccWithPedigreeInfo(canvasSegments, copyNumberModel, pedigreeInfo);
            if (CanvasPedigreeCaller.UseAlleleCountsInformation(canvasSegments, _callerParameters.MinAlleleCountsThreshold, _callerParameters.MinAlleleNumberInSegment) &&
                pedigreeInfo.HasOther())
                AssignMccNoPedigreeInfo(canvasSegments.Where(segment => pedigreeInfo.OtherIds.Contains(segment.SampleId)).ToSampleMap(), copyNumberModel, _genotypes);

        }

        /// <summary>
        /// Calculates maximal likelihood for segments with SNV allele counts given CopyNumber. Updated MajorChromosomeCount.
        /// </summary>   
        private void AssignMccNoPedigreeInfo(ISampleMap<CanvasSegment> canvasSegments,
            ISampleMap<ICopyNumberModel> model, Dictionary<int, List<PhasedGenotype>> genotypes)
        {
            const int diploidCopyNumber = 2;
            foreach (var sampleId in canvasSegments.SampleIds)
            {
                // variant caller does not attempt to call LOH, for DELs CN=MCC
                int copyNumber = canvasSegments[sampleId].CopyNumber;
                if (copyNumber <= diploidCopyNumber)
                {
                    if (copyNumber == diploidCopyNumber)
                        canvasSegments[sampleId].MajorChromosomeCount = null;
                    else
                        canvasSegments[sampleId].MajorChromosomeCount = copyNumber;
                    continue;
                }
                var genotypeset = genotypes[copyNumber];
                int? selectedGtState = null;
                double gqscore = GetGtLogLikelihoodScore(canvasSegments[sampleId].Balleles, genotypeset, ref selectedGtState, model[sampleId]);
                if (selectedGtState.HasValue)
                {
                    canvasSegments[sampleId].MajorChromosomeCount =
                        Math.Max(genotypeset[selectedGtState.Value].CopyNumberA,
                            genotypeset[selectedGtState.Value].CopyNumberB);
                    canvasSegments[sampleId].MajorChromosomeCountScore = gqscore;
                }

            }
        }

        /// <summary>
        /// Calculates maximal likelihood for genotypes given a copy number call. Updated MajorChromosomeCount.
        /// </summary>
        private void AssignMccWithPedigreeInfo(ISampleMap<CanvasSegment> canvasSegments, ISampleMap<ICopyNumberModel> model, PedigreeInfo pedigreeInfo)
        {
            double maximalLogLikelihood = Double.NegativeInfinity;
            int parent1CopyNumber = canvasSegments[pedigreeInfo.ParentsIds.First()].CopyNumber;
            int parent2CopyNumber = canvasSegments[pedigreeInfo.ParentsIds.Last()].CopyNumber;

            foreach (var parent1GtStates in _genotypes[parent1CopyNumber])
            {
                foreach (var parent2GtStates in _genotypes[parent2CopyNumber])
                {
                    var bestChildGtStates = new List<PhasedGenotype>();
                    double currentLogLikelihood = 0;
                    foreach (SampleId child in pedigreeInfo.OffspringIds)
                    {
                        int childCopyNumber = canvasSegments[child].CopyNumber;
                        bool isInheritedCnv = !canvasSegments[child].DqScore.HasValue;
                        double bestLogLikelihood = Double.NegativeInfinity;
                        PhasedGenotype bestGtState = null;
                        bestLogLikelihood = GetProbandLogLikelihood(model[child], childCopyNumber,
                            parent1GtStates, parent2GtStates, isInheritedCnv, canvasSegments[child], bestLogLikelihood, ref bestGtState);
                        bestChildGtStates.Add(bestGtState);
                        currentLogLikelihood += bestLogLikelihood;
                    }
                    currentLogLikelihood += GetCurrentGtLogLikelihood(model[pedigreeInfo.ParentsIds.First()], canvasSegments[pedigreeInfo.ParentsIds.First()], parent1GtStates) +
                                            GetCurrentGtLogLikelihood(model[pedigreeInfo.ParentsIds.Last()], canvasSegments[pedigreeInfo.ParentsIds.Last()], parent2GtStates);

                    currentLogLikelihood = Double.IsNaN(currentLogLikelihood) || Double.IsInfinity(currentLogLikelihood)
                        ? Double.NegativeInfinity
                        : currentLogLikelihood;

                    if (currentLogLikelihood > maximalLogLikelihood)
                    {
                        maximalLogLikelihood = currentLogLikelihood;
                        AssignMcc(canvasSegments[pedigreeInfo.ParentsIds.First()], model[pedigreeInfo.ParentsIds.First()], parent1GtStates, parent1CopyNumber);
                        AssignMcc(canvasSegments[pedigreeInfo.ParentsIds.Last()], model[pedigreeInfo.ParentsIds.Last()], parent2GtStates, parent2CopyNumber);
                        for (int childIndex = 0; childIndex < pedigreeInfo.OffspringIds.Count; childIndex++)
                        {
                            var childId = pedigreeInfo.OffspringIds[childIndex];
                            var bestChildGtState = bestChildGtStates[childIndex];
                            if (bestChildGtState == null) continue;
                            var childSegment = canvasSegments[childId];
                            AssignMcc(childSegment, model[childId], bestChildGtState, childSegment.CopyNumber);
                        }
                    }
                }
            }
        }

        private double GetProbandLogLikelihood(ICopyNumberModel copyNumberModel, int childCopyNumber, PhasedGenotype parent1GtStates, PhasedGenotype parent2GtStates, bool isInheritedCnv, CanvasSegment canvasSegment,
            double bestLogLikelihood, ref PhasedGenotype bestGtState)
        {
            foreach (var childGtState in _genotypes[childCopyNumber])
            {
                double currentChildLogLikelihood;
                if (IsGtPedigreeConsistent(parent1GtStates, childGtState) &&
                    IsGtPedigreeConsistent(parent2GtStates, childGtState)
                    && isInheritedCnv)
                    currentChildLogLikelihood = copyNumberModel.GetGenotypeLogLikelihood(canvasSegment.Balleles, childGtState);
                else
                    continue;
                if (currentChildLogLikelihood > bestLogLikelihood)
                {
                    bestLogLikelihood = currentChildLogLikelihood;
                    bestGtState = childGtState;
                }
            }
            return bestLogLikelihood;
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
            if (copyNumber > diploidCopyNumber)
            {
                canvasSegment.MajorChromosomeCount =
                    Math.Max(gtStates.CopyNumberA, gtStates.CopyNumberB);
                int? selectedGtState = _genotypes[copyNumber].IndexOf(gtStates);
                canvasSegment.MajorChromosomeCountScore = GetGtLogLikelihoodScore(canvasSegment.Balleles, _genotypes[copyNumber], ref selectedGtState, copyNumberModel);
            }
            else
            {
                // variant caller does not attempt to call LOH, for DELs CN=MCC
                if (copyNumber == diploidCopyNumber)
                    canvasSegment.MajorChromosomeCount = null;
                else
                    canvasSegment.MajorChromosomeCount = copyNumber;
                canvasSegment.MajorChromosomeCountScore = null;
            }
        }

        internal static double GetGtLogLikelihoodScore(Balleles gtObservedCounts, List<PhasedGenotype> gtModelCounts, ref int? selectedGtState, ICopyNumberModel copyNumberModel)
        {
            const int maxGQscore = 60;
            var gtLogLikelihoods = Enumerable.Repeat(Double.NegativeInfinity, gtModelCounts.Count).ToList();
            var gtModelCounter = -1;
            foreach (var gtModelCount in gtModelCounts)
            {
                gtModelCounter++;
                // As we don't estimate allele CN but only MCC, focus on upper-triangle 
                if (gtModelCount.CopyNumberA < gtModelCount.CopyNumberB)
                    continue;
                gtLogLikelihoods[gtModelCounter] = copyNumberModel.GetGenotypeLogLikelihood(gtObservedCounts, gtModelCount);
            }
            var maxLogLikelihood = gtLogLikelihoods.Max();
            if (!selectedGtState.HasValue)
                selectedGtState = gtLogLikelihoods.IndexOf(maxLogLikelihood);
            double normalizationConstant = gtLogLikelihoods.Sum(ll => Math.Exp(ll - maxLogLikelihood));
            double gqscore = -10.0 * Math.Log10((normalizationConstant - 1) / normalizationConstant);
            if (Double.IsInfinity(gqscore) | gqscore > maxGQscore)
                gqscore = maxGQscore;
            return Double.IsNaN(gqscore) || Double.IsInfinity(gqscore) ? 0 : gqscore;
        }

        private static double GetCurrentGtLogLikelihood(ICopyNumberModel copyNumberModel, CanvasSegment canvasSegment, PhasedGenotype gtStates)
        {
            return copyNumberModel.GetGenotypeLogLikelihood(canvasSegment.Balleles, gtStates);
        }

        /// <summary>
        /// Estimate joint likelihood and most likely CN assignment within a pedigree using total CN Genotype likelihoods and transition matrix
        /// </summary>
        /// <param name="pedigreeInfo"></param>
        /// <param name="copyNumbersLikelihoods"></param>
        /// <returns></returns>
        private (ISampleMap<Genotype> copyNumbersGenotypes, JointLikelihoods jointLikelihood) GetPedigreeCopyNumbers(PedigreeInfo pedigreeInfo, ISampleMap<Dictionary<Genotype, double>> copyNumbersLikelihoods)
        {
            int nHighestLikelihoodGenotypes = pedigreeInfo != null && pedigreeInfo.OffspringIds.Count >= 2 ? 3 : _callerParameters.MaximumCopyNumber;
            copyNumbersLikelihoods = copyNumbersLikelihoods.SelectValues(l => l.OrderByDescending(kvp => kvp.Value).Take(nHighestLikelihoodGenotypes).ToDictionary());

            var sampleCopyNumbersGenotypes = new SampleMap<Genotype>();
            var jointLikelihood = new JointLikelihoods();
            if (!pedigreeInfo.HasFullPedigree())
                return (sampleCopyNumbersGenotypes, jointLikelihood);
            // parent 1 total CNs and likelihoods
            foreach (var copyNumberParent1 in copyNumbersLikelihoods[pedigreeInfo.ParentsIds.First()])
            {
                // parent 2 total CNs and likelihoods
                foreach (var copyNumberParent2 in copyNumbersLikelihoods[pedigreeInfo.ParentsIds.Last()])
                {
                    // for offspring in addition to querying likelihoods using total CNs, iterate over all possible genotype combination (CopyNumberA/B) for a given
                    // CN and estimate likely transition probabilities using TransitionMatrix
                    foreach (var offspringGtStates in pedigreeInfo.OffspringPhasedGenotypes)
                    {
                        if (!pedigreeInfo.OffspringIds.All(id => copyNumbersLikelihoods[id].ContainsKey(
                            Genotype.Create(Math.Min(offspringGtStates[pedigreeInfo.OffspringIds.IndexOf(id)].PhasedGenotype.CopyNumberA + offspringGtStates[pedigreeInfo.OffspringIds.IndexOf(id)].PhasedGenotype.CopyNumberB,
                                _callerParameters.MaximumCopyNumber - 1)))))
                        {
                            // unavailable total CN
                            continue;
                        }
                        // For a given combination of offspring copy numbers, only the genotypes that result in the maximum likelihood contribute to the final result."
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
                            {pedigreeInfo.ParentsIds.Last(),  copyNumberParent2.Key}
                        };
                        pedigreeInfo.OffspringIds.Zip(totalCopyNumberGenotypes).ForEach(sampleIdGenotypeKvp => genotypesInPedigree.Add(sampleIdGenotypeKvp.Item1, sampleIdGenotypeKvp.Item2));
                        genotypesInPedigree = genotypesInPedigree.OrderBy(pedigreeInfo.AllSampleIds);
                        jointLikelihood.AddJointLikelihood(genotypesInPedigree, currentLikelihood);
                        double currentLogLikelihood = Math.Log(currentLikelihood);
                        if (currentLogLikelihood > jointLikelihood.MaximalLogLikelihood)
                        {
                            jointLikelihood.MaximalLogLikelihood = currentLogLikelihood;
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