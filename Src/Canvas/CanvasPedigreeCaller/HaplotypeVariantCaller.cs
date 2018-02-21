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
            ISampleMap<ICopyNumberModel> copyNumberModel, PedigreeInfo pedigreeInfo)
        {
            var coverageLikelihoods = _copyNumberLikelihoodCalculator.GetCopyNumbersLikelihoods(canvasSegments, samplesInfo, copyNumberModel);
            // if number and properties of SNPs in the segment are above threshold, calculate likelihood from SNPs and merge with
            // coverage likelihood to form merged likelihoods
            int nBalleles = canvasSegments.Values.First().Balleles.Size();
            // If allele information is available (i.e. segment has enough SNPs) merge coverage and alelle likelihood obtained by GetGenotypeLogLikelihoods
            // into singleSampleLikelihoods using JoinLikelihoods function. 
            var singleSampleLikelihoods = CanvasPedigreeCaller.UseAlleleCountsInformation(canvasSegments, 
                _callerParameters.MinAlleleCountsThreshold, _callerParameters.MinAlleleNumberInSegment)
                ? JoinLikelihoods(GetGenotypeLogLikelihoods(canvasSegments, copyNumberModel, _PhasedGenotypes), coverageLikelihoods, nBalleles)
                : coverageLikelihoods;
            // estimate joint likelihood across pedigree samples from singleSampleLikelihoods using either only coverage or coverage + allele counts
            (var pedigreeCopyNumbers, var pedigreeLikelihoods) = GetPedigreeCopyNumbers(pedigreeInfo, singleSampleLikelihoods);

            var nonPedigreeCopyNumbers = CanvasPedigreeCaller.GetNonPedigreeCopyNumbers(canvasSegments, pedigreeInfo, singleSampleLikelihoods);

            var mergedCopyNumbers = nonPedigreeCopyNumbers.Concat(pedigreeCopyNumbers).OrderBy(canvasSegments.SampleIds);

            AssignCNandScores(canvasSegments, samplesInfo, pedigreeInfo, singleSampleLikelihoods,
                pedigreeLikelihoods, mergedCopyNumbers);
        }

        private static ISampleMap<Dictionary<PhasedGenotype, double>> GetGenotypeLogLikelihoods(ISampleMap<CanvasSegment> canvasSegments,
            ISampleMap<ICopyNumberModel> copyNumberModel, List<PhasedGenotype> genotypes)
        {
            var REF = new PhasedGenotype(1, 1);
            var loh = new List<PhasedGenotype> { new PhasedGenotype(0, 2), new PhasedGenotype(2, 0) };

            var singleSampleLikelihoods = new SampleMap<Dictionary<PhasedGenotype, double>>();
            foreach (var sampleId in canvasSegments.SampleIds)
            {
                   var logLikelihoods = genotypes.Select(genotype => (genotype, copyNumberModel[sampleId].
                GetGenotypeLogLikelihood(canvasSegments[sampleId].Balleles, genotype))).ToDictionary(kvp => kvp.Item1, kvp => kvp.Item2);
                if (logLikelihoods[REF] >= Math.Max(logLikelihoods[loh.First()], logLikelihoods[loh.Last()]))
                    logLikelihoods[loh.First()] = logLikelihoods[loh.Last()] = logLikelihoods.Values.Where(ll => ll > Double.NegativeInfinity).Min();
                singleSampleLikelihoods.Add(sampleId, logLikelihoods);
            }
            return singleSampleLikelihoods;
        }

        /// <summary>
        /// Merge likelihoods separately derived from coverage or allele counts (from SNPs) data
        /// </summary>
        /// <param name="genotypeLogLikelihoods"></param>
        /// <param name="copyNumberLikelihoods"></param>
        /// <returns></returns>
        private ISampleMap<Dictionary<Genotype, double>> JoinLikelihoods(ISampleMap<Dictionary<PhasedGenotype, double>> genotypeLogLikelihoods,
            ISampleMap<Dictionary<Genotype, double>> copyNumberLikelihoods, int nBalleles)
        {
            double minLogLikelihood = double.MinValue;
            var jointSampleLogLikelihoods = new SampleMap<Dictionary<Genotype, double>>();

            foreach (var sampleId in genotypeLogLikelihoods.SampleIds)
            {
                var jointLogLikelihoods = new Dictionary<Genotype, double>();
                // since both likelihoods use negative binomial distribution, area under the curve should be 
                // similar except for all b-alleles vs median point estimate used for depth, so correct for this
                // by using number of alleles
                foreach (var genotypeLikelihood in genotypeLogLikelihoods[sampleId])
                {
                    int totalCopyNumber = genotypeLikelihood.Key.CopyNumberA + genotypeLikelihood.Key.CopyNumberB;
                    double jointLogLikelihood = genotypeLikelihood.Value / nBalleles +
                        Math.Max(minLogLikelihood, Math.Log(copyNumberLikelihoods[sampleId][Genotype.Create(totalCopyNumber)]));
                    jointLogLikelihoods[Genotype.Create(genotypeLikelihood.Key)] = jointLogLikelihood;
                }
                jointSampleLogLikelihoods.Add(sampleId, jointLogLikelihoods);
            }
            return jointSampleLogLikelihoods;
        }

        /// <summary>
        /// Evaluate joint likelihood of all genotype combinations across samples. 
        /// Return joint likelihood object and the copy number states with the highest likelihood 
        /// </summary>
        private (SampleMap<Genotype> copyNumbersGenotypes, JointLikelihoods jointLikelihood) GetPedigreeCopyNumbers(PedigreeInfo pedigreeInfo, ISampleMap<Dictionary<Genotype, double>> singleSampleLogLikelihoods)
        {
            int nHighestLikelihoodGenotypes = pedigreeInfo != null && pedigreeInfo.OffspringIds.Count >= 2 ? 3 : _callerParameters.MaximumCopyNumber;
            singleSampleLogLikelihoods = singleSampleLogLikelihoods.SelectValues(l => l.OrderByDescending(kvp => kvp.Value).Take(nHighestLikelihoodGenotypes).ToDictionary());

            var sampleCopyNumbersGenotypes = new SampleMap<Genotype>();
            var jointLikelihood = new JointLikelihoods();
            if (!pedigreeInfo.HasFullPedigree())
                return (sampleCopyNumbersGenotypes, jointLikelihood);
            var usePhasedGenotypes = singleSampleLogLikelihoods.Values.First().Keys.First().PhasedGenotype != null;
            double maxOffspringLogLikelihood = pedigreeInfo.OffspringIds.Select(key => singleSampleLogLikelihoods[key].Select(kvp => kvp.Value).Max()).Aggregate(1.0, (acc, val) => acc * val);

            foreach (var parent1GtStates in singleSampleLogLikelihoods[pedigreeInfo.ParentsIds.First()])
            {
                foreach (var parent2GtStates in singleSampleLogLikelihoods[pedigreeInfo.ParentsIds.Last()])
                {
                    foreach (var genotypes in usePhasedGenotypes ? pedigreeInfo.OffspringPhasedGenotypes : pedigreeInfo.OffspringTotalCopyNumberGenotypes)
                    {

                        double currentLogLikelihood = parent1GtStates.Value + parent2GtStates.Value;
                        if (currentLogLikelihood + maxOffspringLogLikelihood <= jointLikelihood.MaximalLogLikelihood)
                        {
                            continue;
                        }
                        if (!pedigreeInfo.OffspringIds.All(id => singleSampleLogLikelihoods[id].ContainsKey(genotypes[pedigreeInfo.OffspringIds.IndexOf(id)])))
                        {
                            continue;
                        }

                        var offspringGtStates = new List<Genotype>();
                        for (int index = 0; index < pedigreeInfo.OffspringIds.Count; index++)
                        {
                            var offspringId = pedigreeInfo.OffspringIds[index];
                            double tmpLikelihood = singleSampleLogLikelihoods[offspringId][genotypes[index]];
                            offspringGtStates.Add(genotypes[index]);
                            currentLogLikelihood += tmpLikelihood;
                            currentLogLikelihood += EstimateTransmissionProbability(parent1GtStates, parent2GtStates,
                                new KeyValuePair<Genotype, double>(genotypes[index], tmpLikelihood), _callerParameters.DeNovoRate, pedigreeInfo);
                        }
                        currentLogLikelihood = Double.IsNaN(currentLogLikelihood) || Double.IsInfinity(currentLogLikelihood) ? Double.MinValue : currentLogLikelihood;
                        var genotypesInPedigree = new SampleMap<Genotype>
                        {
                            {pedigreeInfo.ParentsIds.First(), parent1GtStates.Key},
                            {pedigreeInfo.ParentsIds.Last(), parent2GtStates.Key}
                        };
                        pedigreeInfo.OffspringIds.Zip(offspringGtStates).ForEach(sampleIdGenotypeKvp => genotypesInPedigree.Add(sampleIdGenotypeKvp.Item1, sampleIdGenotypeKvp.Item2));
                        // return to SampleId ordering 
                        var orderedGenotypesInPedigree = genotypesInPedigree.OrderBy(x => singleSampleLogLikelihoods.SampleIds.ToList().IndexOf(x.Key)).ToSampleMap();
                        // convert to likelihood
                        jointLikelihood.AddJointLikelihood(orderedGenotypesInPedigree, Double.IsNaN(Math.Exp(currentLogLikelihood)) ? 0 : Math.Exp(currentLogLikelihood));

                        if (currentLogLikelihood > jointLikelihood.MaximalLogLikelihood)
                        {
                            jointLikelihood.MaximalLogLikelihood = currentLogLikelihood;
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
                    GetSingleSampleQualityScore(singleSampleLikelihoods[sampleId], copyNumbers[sampleId], pedigreeMembersInfo[sampleId].GetPloidy(canvasSegments[sampleId]));
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

        private double GetSingleSampleQualityScore(Dictionary<Genotype, double> copyNumbersLikelihoods, Genotype selectedGenotype, int ploidy)
        {
            var altCn = copyNumbersLikelihoods.Keys.Where(gt => gt.TotalCopyNumber == selectedGenotype.TotalCopyNumber);
            var refCn = selectedGenotype.TotalCopyNumber != ploidy ? copyNumbersLikelihoods.Keys.Where(gt => gt.TotalCopyNumber == ploidy) :
                copyNumbersLikelihoods.Keys.Where(gt => gt.TotalCopyNumber != selectedGenotype.TotalCopyNumber);
            var allCn = altCn.Concat(refCn);

            double maxLogLikelihood = copyNumbersLikelihoods.Max(ll => ll.Value);
            double normalizationConstant = copyNumbersLikelihoods.SelectKeys(allCn).Sum(ll => Math.Exp(ll.Value - maxLogLikelihood));
            double altLikelihood = copyNumbersLikelihoods.SelectKeys(altCn).Sum(ll => Math.Exp(ll.Value - maxLogLikelihood));
            double qscore = -10.0 * Math.Log10((normalizationConstant - altLikelihood) / normalizationConstant);
            if (Double.IsInfinity(qscore) | qscore > _callerParameters.MaxQscore)
                qscore = _callerParameters.MaxQscore;
            return qscore;
        }

        /// <summary>
        /// Assess likelihood of a de-novo variant for copyNumberGenotypes configuration with a Mendelian conflict 
        /// </summary>
        /// <param name="jointLikelihoods"></param>
        /// <param name="probandId"></param>
        /// <param name="copyNumberGenotypes"></param>
        /// <returns></returns>
        private double GetConditionalDeNovoQualityScore(JointLikelihoods jointLikelihoods, SampleId probandId, ISampleMap<Genotype> copyNumberGenotypes)
        {
            const double q60 = 0.000001;
            // likelihood of proband having copyNumberGenotypes in a pedigree
            double marginalAltLikelihood = jointLikelihoods.GetMarginalLikelihood(new KeyValuePair<SampleId, Genotype>(probandId, copyNumberGenotypes[probandId]));
            // likelihood of proband not having copyNumberGenotypes in a pedigree
            double marginalRefLikelihood = jointLikelihoods.GetMarginalNonAltLikelihood(new KeyValuePair<SampleId, Genotype>(probandId, copyNumberGenotypes[probandId]));
            double normalization = marginalAltLikelihood + marginalRefLikelihood;
            // likelihood of proband genotype != ALT given "copyNumberGenotypes" configuration in pedigree with Mendelian conflict 
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