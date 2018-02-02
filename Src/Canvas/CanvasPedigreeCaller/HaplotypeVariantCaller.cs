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
    internal class HaplotypeVariantCaller : IVariantCaller
    {
        private readonly CopyNumberLikelihoodCalculator _copyNumberLikelihoodCalculator;
        private readonly List<PhasedGenotype> _PhasedGenotypes;
        private readonly int _qualityFilterThreshold;
        private readonly PedigreeCallerParameters _callerParameters;

        public HaplotypeVariantCaller(CopyNumberLikelihoodCalculator copyNumberLikelihoodCalculator, PedigreeCallerParameters callerParameters, int qualityFilterThreshold)
        {
            _copyNumberLikelihoodCalculator = copyNumberLikelihoodCalculator;
            _PhasedGenotypes = GenerateGenotypeCombinations(callerParameters.MaximumCopyNumber);
            _callerParameters = callerParameters;
            _qualityFilterThreshold = qualityFilterThreshold;
        }

        public void CallVariant(ISampleMap<CanvasSegment> canvasSegments, ISampleMap<SampleMetrics> samplesInfo,
            ISampleMap<ICopyNumberModel> copyNumberModel,
            PedigreeInfo pedigreeInfo)
        {
            var coverageLikelihoods = _copyNumberLikelihoodCalculator.GetCopyNumbersLikelihoods(canvasSegments, samplesInfo, copyNumberModel);
            // if number and properties of SNPs in the segment are above threshold, calculate likelihood from SNPs and merge with
            // coverage likelihood to form merged likelihoods
            var singleSampleLikelihoods = UseAlleleCountsInformation(canvasSegments)
                ? JoinLikelihoods(GetGenotypeLikelihoods(canvasSegments, copyNumberModel, _PhasedGenotypes), coverageLikelihoods)
                : coverageLikelihoods;

            (var pedigreeCopyNumbers, var pedigreeLikelihoods) = GetPedigreeCopyNumbers(pedigreeInfo, singleSampleLikelihoods);

            var nonPedigreeCopyNumbers = CanvasPedigreeCaller.GetNonPedigreeCopyNumbers(canvasSegments, pedigreeInfo, singleSampleLikelihoods);

            var mergedCopyNumbers = nonPedigreeCopyNumbers.Concat(pedigreeCopyNumbers).OrderBy(canvasSegments.SampleIds);

            AssignCNandScores(canvasSegments, samplesInfo, pedigreeInfo, singleSampleLikelihoods,
                pedigreeLikelihoods, mergedCopyNumbers);
        }


        /// <summary>
        /// Derives metrics from b-allele counts within each segment and determines whereas to use them for calculating MCC
        /// </summary>
        /// <param name="canvasSegments"></param>
        /// <returns></returns>
        private bool UseAlleleCountsInformation(ISampleMap<CanvasSegment> canvasSegments)
        {
            var alleles = canvasSegments.Values.Select(segments => segments.Balleles?.TotalCoverage);
            var alleleCounts = alleles.Select(allele => allele?.Count ?? 0).ToList();
            bool lowAlleleCounts = alleleCounts.Select(x => x < _callerParameters.DefaultReadCountsThreshold).Any(c => c);
            // var coverageCounts = canvasSegments.Values.Select(segments => segments.MedianCount).ToList();
            // double alleleDensity = canvasSegments.Values.First().Length / Math.Max(alleleCounts.Average(), 1.0);
            // bool useCnLikelihood = lowAlleleCounts || alleleDensity < _callerParameters.DefaultAlleleDensityThreshold || alleleCounts.Any(x => x > _callerParameters.DefaultPerSegmentAlleleMaxCounts) || coverageCounts.Any(coverage => coverage < medianCoverageThreshold); 
            // for now only use lowAlleleCounts metric
            return !lowAlleleCounts;
        }

        private ISampleMap<Dictionary<Genotype, double>> GetCopyNumbersLikelihoods(ISampleMap<CanvasSegment> canvasSegments, ISampleMap<SampleMetrics> samplesInfo,
            ISampleMap<CopyNumberModel> copyNumberModel, List<Genotype> genotypes)
        {
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

        private static ISampleMap<Dictionary<PhasedGenotype, double>> GetGenotypeLikelihoods(ISampleMap<CanvasSegment> canvasSegments,
            ISampleMap<ICopyNumberModel> copyNumberModel, List<PhasedGenotype> genotypes)
        {
            var REF = new PhasedGenotype(1, 1);
            var loh = new List<PhasedGenotype> { new PhasedGenotype(0, 2), new PhasedGenotype(2, 0) };

            var singleSampleLikelihoods = new SampleMap<Dictionary<PhasedGenotype, double>>();
            foreach (var sampleId in canvasSegments.SampleIds)
            {
                var likelihoods = genotypes.Select(genotype => (genotype, copyNumberModel[sampleId].GetGenotypeLikelihood(canvasSegments[sampleId].Balleles, genotype))).ToDictionary(kvp => kvp.Item1, kvp => kvp.Item2);
                if (likelihoods[REF] >= Math.Max(likelihoods[loh.First()], likelihoods[loh.Last()]))
                    likelihoods[loh.First()] = likelihoods[loh.Last()] = likelihoods.Values.Where(ll => ll > 0).Min();
                else if (likelihoods[loh.First()] > Math.Max(likelihoods[REF], likelihoods[loh.Last()]))
                    likelihoods[REF] = likelihoods[loh.Last()] = likelihoods.Values.Where(ll => ll > 0).Min();
                else if(likelihoods[loh.Last()] > Math.Max(likelihoods[loh.First()], likelihoods[REF]))
                    likelihoods[REF] = likelihoods[loh.First()] = likelihoods.Values.Where(ll => ll > 0).Min();
                singleSampleLikelihoods.Add(sampleId, likelihoods);
            }
            return singleSampleLikelihoods;
        }

        /// <summary>
        /// Merge likelihoods separately derived from coverage or allele counts (from SNPs) data
        /// </summary>
        /// <param name="genotypeLikelihoods"></param>
        /// <param name="copyNumberLikelihoods"></param>
        /// <returns></returns>
        private ISampleMap<Dictionary<Genotype, double>> JoinLikelihoods(ISampleMap<Dictionary<PhasedGenotype, double>> genotypeLikelihoods,
            ISampleMap<Dictionary<Genotype, double>> copyNumberLikelihoods)
        {
            var jointSampleLikelihoods = new SampleMap<Dictionary<Genotype, double>>();
            foreach (var sampleId in genotypeLikelihoods.SampleIds)
            {
                var jointLikelihoods = new Dictionary<Genotype, double>();
                double areaGt = genotypeLikelihoods[sampleId].Values.Sum();
                double areaCn = copyNumberLikelihoods[sampleId].Values.Sum();
                foreach (var genotypeLikelihood in genotypeLikelihoods[sampleId])
                {
                    double scaler = areaGt > areaCn ? areaCn / areaGt : areaGt / areaCn;
                    int totalCopyNumber = genotypeLikelihood.Key.CopyNumberA + genotypeLikelihood.Key.CopyNumberB;
                    double jointLikelihood = genotypeLikelihood.Value * scaler *
                                             copyNumberLikelihoods[sampleId][Genotype.Create(totalCopyNumber)];
                    jointLikelihoods[Genotype.Create(genotypeLikelihood.Key)] = jointLikelihood;
                }
                jointSampleLikelihoods.Add(sampleId, jointLikelihoods);
            }
            return jointSampleLikelihoods;
        }

        /// <summary>
        /// Evaluate joint likelihood of all genotype combinations across samples. 
        /// Return joint likelihood object and the copy number states with the highest likelihood 
        /// </summary>
        private (SampleMap<Genotype> copyNumbersGenotypes, JointLikelihoods jointLikelihood) GetPedigreeCopyNumbers(PedigreeInfo pedigreeInfo, ISampleMap<Dictionary<Genotype, double>> singleSampleLikelihoods)
        {
            int nHighestLikelihoodGenotypes = pedigreeInfo != null && pedigreeInfo.OffspringIds.Count >= 2 ? 3 : _callerParameters.MaximumCopyNumber;
            singleSampleLikelihoods = singleSampleLikelihoods.SelectValues(l => l.OrderByDescending(kvp => kvp.Value).Take(nHighestLikelihoodGenotypes).ToDictionary());

            var sampleCopyNumbersGenotypes = new SampleMap<Genotype>();
            var jointLikelihood = new JointLikelihoods();
            if (!pedigreeInfo.HasFullPedigree())
                return (sampleCopyNumbersGenotypes, jointLikelihood);
            var usePhasedGenotypes = singleSampleLikelihoods.Values.First().Keys.First().PhasedGenotype != null;
            double maxOffspringLikelihood = pedigreeInfo.OffspringIds.Select(key => singleSampleLikelihoods[key].Select(kvp => kvp.Value).Max()).Aggregate(1.0, (acc, val) => acc * val);

            foreach (var parent1GtStates in singleSampleLikelihoods[pedigreeInfo.ParentsIds.First()])
            {
                foreach (var parent2GtStates in singleSampleLikelihoods[pedigreeInfo.ParentsIds.Last()])
                {
                    foreach (var genotypes in usePhasedGenotypes ? pedigreeInfo.OffspringPhasedGenotypes : pedigreeInfo.OffspringTotalCopyNumberGenotypes)
                    {

                        double currentLikelihood = parent1GtStates.Value * parent2GtStates.Value;
                        if (currentLikelihood * maxOffspringLikelihood <= jointLikelihood.MaximalLikelihood)
                        {
                            continue;
                        }
                        if (!pedigreeInfo.OffspringIds.All(id => singleSampleLikelihoods[id].ContainsKey(genotypes[pedigreeInfo.OffspringIds.IndexOf(id)])))
                        {
                            continue;
                        }

                        var offspringGtStates = new List<Genotype>();
                        for (int index = 0; index < pedigreeInfo.OffspringIds.Count; index++)
                        {
                            var offspringId = pedigreeInfo.OffspringIds[index];
                            double tmpLikelihood = singleSampleLikelihoods[offspringId][genotypes[index]];
                            offspringGtStates.Add(genotypes[index]);
                            currentLikelihood *= tmpLikelihood;
                            currentLikelihood *= EstimateTransmissionProbability(parent1GtStates, parent2GtStates,
                                new KeyValuePair<Genotype, double>(genotypes[index], tmpLikelihood), _callerParameters.DeNovoRate, pedigreeInfo);
                        }
                        currentLikelihood = Double.IsNaN(currentLikelihood) || Double.IsInfinity(currentLikelihood) ? 0 : currentLikelihood;
                        var genotypesInPedigree = new SampleMap<Genotype>
                        {
                            {pedigreeInfo.ParentsIds.First(), parent1GtStates.Key},
                            {pedigreeInfo.ParentsIds.Last(), parent2GtStates.Key}
                        };
                        pedigreeInfo.OffspringIds.Zip(offspringGtStates).ForEach(sampleIdGenotypeKvp => genotypesInPedigree.Add(sampleIdGenotypeKvp.Item1, sampleIdGenotypeKvp.Item2));
                        // return to SampleId ordering 
                        var orderedGenotypesInPedigree = genotypesInPedigree.OrderBy(x => singleSampleLikelihoods.SampleIds.ToList().IndexOf(x.Key)).ToSampleMap();
                        jointLikelihood.AddJointLikelihood(orderedGenotypesInPedigree, currentLikelihood);

                        if (currentLikelihood > jointLikelihood.MaximalLikelihood)
                        {
                            jointLikelihood.MaximalLikelihood = currentLikelihood;
                            sampleCopyNumbersGenotypes = orderedGenotypesInPedigree;
                        }
                    }
                }
            }
            if (sampleCopyNumbersGenotypes.Empty())
                throw new IlluminaException("Maximal likelihood was not found");
            return (sampleCopyNumbersGenotypes, jointLikelihood);
        }

        /// <summary>
        /// Estimate Transmission probability for parental copy number genotypes.
        /// Uses de novo rate when genotypes can be evaluated (segment with SNP).
        /// </summary>
        /// <param name="parent1GtStates"></param>
        /// <param name="parent2GtStates"></param>
        /// <param name="offspringGtState"></param>
        /// <param name="deNovoRate"></param>
        /// <param name="pedigreeInfo"></param>
        /// <returns></returns>
        private static double EstimateTransmissionProbability(KeyValuePair<Genotype, double> parent1GtStates, KeyValuePair<Genotype, double> parent2GtStates, KeyValuePair<Genotype, double> offspringGtState, double deNovoRate, PedigreeInfo pedigreeInfo)
        {
            if (parent1GtStates.Key.HasAlleleCopyNumbers && parent2GtStates.Key.HasAlleleCopyNumbers)
                return (offspringGtState.Key.PhasedGenotype.ContainsSharedAlleleA(parent1GtStates.Key.PhasedGenotype) ||
                       offspringGtState.Key.PhasedGenotype.ContainsSharedAlleleA(parent2GtStates.Key.PhasedGenotype)) &&
                       (offspringGtState.Key.PhasedGenotype.ContainsSharedAlleleB(parent1GtStates.Key.PhasedGenotype) ||
                        offspringGtState.Key.PhasedGenotype.ContainsSharedAlleleB(parent2GtStates.Key.PhasedGenotype))
                    ? 1.0
                    : deNovoRate;
            return pedigreeInfo.TransitionMatrix[parent1GtStates.Key.TotalCopyNumber][
                       offspringGtState.Key.TotalCopyNumber] *
                   pedigreeInfo.TransitionMatrix[parent2GtStates.Key.TotalCopyNumber][
                       offspringGtState.Key.TotalCopyNumber];
        }

        private void AssignCNandScores(ISampleMap<CanvasSegment> canvasSegments,
            ISampleMap<SampleMetrics> pedigreeMembersInfo,
            PedigreeInfo pedigreeInfo, ISampleMap<Dictionary<Genotype, double>> singleSampleLikelihoods,
            JointLikelihoods jointLikelihoods, ISampleMap<Genotype> copyNumbers)
        {
            foreach (var sampleId in canvasSegments.SampleIds)
            {
                canvasSegments[sampleId].QScore =
                    GetSingleSampleQualityScore(singleSampleLikelihoods[sampleId], copyNumbers[sampleId]);
                canvasSegments[sampleId].CopyNumber = copyNumbers[sampleId].TotalCopyNumber;
                if (canvasSegments[sampleId].QScore < _qualityFilterThreshold)
                    canvasSegments[sampleId].Filter = CanvasFilter.Create(new[] {$"q{_qualityFilterThreshold}"});
                if (copyNumbers[sampleId].PhasedGenotype != null)
                    canvasSegments[sampleId].MajorChromosomeCount =
                        Math.Max(copyNumbers[sampleId].PhasedGenotype.CopyNumberA,
                            copyNumbers[sampleId].PhasedGenotype.CopyNumberB);
            }
            if (pedigreeInfo.HasFullPedigree())
            {
                var pedigreeMembers = pedigreeInfo.ParentsIds.Concat(pedigreeInfo.OffspringIds).ToList();
                var pedigreeMemberCopyNumbers = copyNumbers.WhereSampleIds(sampleId => pedigreeMembers.Contains(sampleId));
                SetDenovoQualityScores(canvasSegments, pedigreeMembersInfo, pedigreeInfo.ParentsIds, pedigreeInfo.OffspringIds, jointLikelihoods, pedigreeMemberCopyNumbers);
            }
    }

        private void SetDenovoQualityScores(ISampleMap<CanvasSegment> canvasSegments, ISampleMap<SampleMetrics> samplesInfo, List<SampleId> parentIDs, List<SampleId> offspringIDs,
            JointLikelihoods jointLikelihoods, ISampleMap<Genotype> copyNumbers)
        {

            foreach (var probandId in offspringIDs)
            {
                // targeted proband is REF
                if (IsReferenceVariant(canvasSegments, samplesInfo, probandId))
                    continue;
                // common variant
                if (IsCommonCnv(copyNumbers, parentIDs, probandId))
                    continue;
                // other offsprings are ALT
                if (!offspringIDs.Except(probandId.ToEnumerable()).All(id => IsReferenceVariant(canvasSegments, samplesInfo, id)))
                    continue;
                // not all q-scores are above the threshold
                if (parentIDs.Concat(probandId).Any(id => !IsPassVariant(canvasSegments, id)))
                    continue;

                double deNovoQualityScore = GetConditionalDeNovoQualityScore(jointLikelihoods, probandId, copyNumbers);
                if (Double.IsInfinity(deNovoQualityScore) | deNovoQualityScore > _callerParameters.MaxQscore)
                    deNovoQualityScore = _callerParameters.MaxQscore;
                canvasSegments[probandId].DqScore = deNovoQualityScore;
            }
        }

        private bool IsPassVariant(ISampleMap<CanvasSegment> canvasSegments, SampleId sampleId)
        {
            return canvasSegments[sampleId].QScore > _qualityFilterThreshold;
        }

        public static bool IsCommonCnv(ISampleMap<Genotype> copyNumberGenotypes, List<SampleId> parentIDs, SampleId probandId)
        {
            var proband = copyNumberGenotypes[probandId];
            var parent1 = copyNumberGenotypes[parentIDs.First()];
            var parent2 = copyNumberGenotypes[parentIDs.Last()];

            if (proband.PhasedGenotype == null)
                return proband.ContainsSharedAlleles(parent1) || proband.ContainsSharedAlleles(parent2) ||
                       proband.TotalCopyNumber > 3 && parent1.TotalCopyNumber > 3 && parent2.TotalCopyNumber > 3 ||
                       proband.TotalCopyNumber == 0 && parent1.TotalCopyNumber == 1 && parent2.TotalCopyNumber == 1;

            return (proband.PhasedGenotype.ContainsSharedAlleleA(parent1.PhasedGenotype) ||
                    proband.PhasedGenotype.ContainsSharedAlleleA(parent2.PhasedGenotype)) &&
                   (proband.PhasedGenotype.ContainsSharedAlleleB(parent1.PhasedGenotype) ||
                    proband.PhasedGenotype.ContainsSharedAlleleB(parent2.PhasedGenotype));
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

        public static int AggregateVariantCoverage(ref List<CanvasSegment> segments)
        {
            var variantCoverage = segments.SelectMany(segment => segment.Balleles.TotalCoverage).ToList();
            return variantCoverage.Any() ? Utilities.Median(variantCoverage) : 0;
        }

        private double GetSingleSampleQualityScore(Dictionary<Genotype, double> copyNumbersLikelihoods, Genotype selectedGenotype)
        {
            double normalizationConstant = copyNumbersLikelihoods.Select(ll => ll.Value).Sum();
            var selectedGenotypes = copyNumbersLikelihoods.Keys.Where(gt => gt.TotalCopyNumber == selectedGenotype.TotalCopyNumber);
            var altLikelihood = copyNumbersLikelihoods.SelectKeys(selectedGenotypes).Select(ll => ll.Value).Sum();
            double qscore = -10.0 * Math.Log10((normalizationConstant - altLikelihood) / normalizationConstant);
            if (Double.IsInfinity(qscore) | qscore > _callerParameters.MaxQscore)
                qscore = _callerParameters.MaxQscore;
            return qscore;
        }

        private double GetConditionalDeNovoQualityScore(JointLikelihoods jointLikelihoods, SampleId probandId, ISampleMap<Genotype> copyNumberGenotypes)
        {
            const double q60 = 0.000001;
            double marginalAltLikelihood = jointLikelihoods.GetMarginalLikelihood(new KeyValuePair<SampleId, Genotype>(probandId, copyNumberGenotypes[probandId]));
            double marginalRefLikelihood = jointLikelihoods.GetMarginalNonAltLikelihood(new KeyValuePair<SampleId, Genotype>(probandId, copyNumberGenotypes[probandId]));
            double normalization = marginalAltLikelihood + marginalRefLikelihood;
            double denovoProbability = jointLikelihoods.GetJointLikelihood(copyNumberGenotypes) / marginalAltLikelihood * marginalRefLikelihood / normalization;
            return -10.0 * Math.Log10(Math.Max(denovoProbability, q60));
        }


        /// <summary>
        /// Generate all possible copy number genotype combinations with the maximal number of alleles per segment set to maxAlleleNumber.
        /// </summary>
        /// <param name="numberOfCnStates"></param>
        /// <returns> </returns>
        public static List<PhasedGenotype> GenerateGenotypeCombinations(int numberOfCnStates)
        {
            var genotypes = new List<PhasedGenotype>();
            for (int cn = 0; cn < numberOfCnStates; cn++)
            {
                for (int gt = 0; gt <= cn; gt++)
                {
                    genotypes.Add(new PhasedGenotype(gt, cn - gt));
                }
            }
            return genotypes;
        }
    }
}