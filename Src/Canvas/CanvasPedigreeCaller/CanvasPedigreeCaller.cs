using CanvasCommon;
using Combinatorics.Collections;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Isas.Framework.DataTypes;
using Isas.Framework.Logging;
using Isas.SequencingFiles;
using Isas.SequencingFiles.Vcf;
using MathNet.Numerics.Distributions;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using Genotype = CanvasCommon.Genotype;


namespace CanvasPedigreeCaller
{

    internal class SegmentIndexRange
    {
        public int Start { get; }
        public int End { get; }

        public SegmentIndexRange(int start, int end)
        {
            Start = start;
            End = end;
        }
    }


    class CanvasPedigreeCaller
    {
        #region Members

        public int QualityFilterThreshold { get; set; } = 7;
        public int DeNovoQualityFilterThreshold { get; set; } = 20;
        public PedigreeCallerParameters CallerParameters { get; set; }
        protected double MedianCoverageThreshold = 4;
        private readonly ILogger _logger;

        public CanvasPedigreeCaller(ILogger logger)
        {
            _logger = logger;
        }

        #endregion

        internal int CallVariants(List<string> variantFrequencyFiles, List<string> segmentFiles,
            string outVcfFile, string ploidyBedPath, string referenceFolder, List<string> sampleNames, string commonCNVsbedPath, string pedigreeFile)
        {
            // load files
            // initialize data structures and classes
            var fileCounter = 0;
            var samplesInfo = new SampleList<SamplesInfo>();
            var sampleSegments = new SampleList<Segments>();
            var copyNumberModels = new SampleList<CopyNumberModel>();
            var variantFrequencyFilesSampleList = new SampleList<string>();
            var kinships = new SampleList<Kinship>();
            if (pedigreeFile != null)
            {
                kinships = ReadPedigreeFile(pedigreeFile);
                // In Kinship enum Proband gets the highest int value 
                var newSampleNames = kinships.OrderByDescending(x => x.Value).Select(x => x.Key.ToString()).ToList();
                var remapIndices = newSampleNames.Select(newname => sampleNames.FindIndex(name => name == newname))
                    .ToList();
                segmentFiles = remapIndices.Select(index => segmentFiles[index]).ToList();
                variantFrequencyFiles = remapIndices.Select(index => variantFrequencyFiles[index]).ToList();
                sampleNames = newSampleNames;               
            }
            else
            {
                sampleNames.ForEach(name => kinships.Add(new SampleId(name), Kinship.Other));
            }

            foreach (string sampleName in sampleNames)
            {
                var sampleId = new SampleId(sampleName);
                var segment = Segments.ReadSegments(_logger, new FileLocation(segmentFiles[fileCounter]));
                segment.AddAlleles(CanvasIO.ReadFrequenciesWrapper(_logger, new FileLocation(variantFrequencyFiles[fileCounter]), segment.IntervalsByChromosome));
                sampleSegments.Add(sampleId, segment);
                var sampleInfo = SamplesInfo.GetSampleInfo(segment, ploidyBedPath, CallerParameters.NumberOfTrimmedBins, sampleId);
                var copyNumberModel = new CopyNumberModel(CallerParameters.MaximumCopyNumber, sampleInfo.MeanMafCoverage, sampleInfo.MeanCoverage, sampleInfo.MaxCoverage);
                samplesInfo.Add(sampleId, sampleInfo);
                copyNumberModels.Add(sampleId, copyNumberModel);
                variantFrequencyFilesSampleList.Add(sampleId, variantFrequencyFiles[fileCounter]);
                fileCounter++;
            }
            var segmentSetsFromCommonCnvs = CreateSegmentSetsFromCommonCnvs(variantFrequencyFilesSampleList,
                CallerParameters.DefaultReadCountsThreshold, commonCNVsbedPath, sampleSegments);

            var segmentsForVariantCalling = GetHighestLikelihoodSegments(segmentSetsFromCommonCnvs, samplesInfo, copyNumberModels);
            var genotypes = GenerateGenotypeCombinations(CallerParameters.MaximumCopyNumber);
            PedigreeInfo pedigreeInfo = null;
            if (kinships.SampleData.Any(kin => kin == Kinship.Proband))
                pedigreeInfo = new PedigreeInfo(kinships, CallerParameters);
            Parallel.ForEach(
                segmentsForVariantCalling,
                new ParallelOptions
                {
                    MaxDegreeOfParallelism = Math.Min(Environment.ProcessorCount, CallerParameters.MaxCoreNumber)
                },
                segments => CallVariant(segments, samplesInfo, copyNumberModels, pedigreeInfo, genotypes)
            );
            var variantCalledSegments = new SampleList<List<CanvasSegment>>();
            foreach (var key in samplesInfo.SampleIds)
                variantCalledSegments.Add(key, segmentsForVariantCalling.Select(segment => segment[key]).ToList());

            var mergedVariantCalledSegments = MergeSegments(variantCalledSegments, CallerParameters.MinimumCallSize);
            var outputFolder = new FileLocation(outVcfFile).Directory;
            foreach (var sampleId in samplesInfo.SampleIds)
            {
                var coverageOutputPath = SingleSampleCallset.GetCoverageAndVariantFrequencyOutput(outputFolder,
                    sampleId.ToString());
                CanvasSegment.WriteCoveragePlotData(mergedVariantCalledSegments[sampleId], samplesInfo[sampleId].MeanCoverage,
                    samplesInfo[sampleId].Ploidy, coverageOutputPath, referenceFolder);
            }
            bool isPedigreeInfoSupplied = pedigreeInfo != null;
            var denovoQualityThreshold = isPedigreeInfoSupplied ? (int?)DeNovoQualityFilterThreshold : null;
            var ploidies = samplesInfo.Select(info => info.Value.Ploidy).ToList();
            var diploidCoverage = samplesInfo.Select(info => info.Value.MeanCoverage).ToList();
            var names = samplesInfo.SampleIds.Select(id => id.ToString()).ToList();
            CanvasSegmentWriter.WriteMultiSampleSegments(outVcfFile, mergedVariantCalledSegments, diploidCoverage, referenceFolder, names,
                null, ploidies, QualityFilterThreshold, isPedigreeInfoSupplied, denovoQualityThreshold);

            outputFolder = new FileLocation(outVcfFile).Directory;
            foreach (var sampleId in samplesInfo.SampleIds)
            {
                var outputVcfPath = SingleSampleCallset.GetSingleSamplePedigreeVcfOutput(outputFolder, sampleId.ToString());
                CanvasSegmentWriter.WriteSegments(outputVcfPath.FullName, mergedVariantCalledSegments[sampleId],
                    samplesInfo[sampleId].MeanCoverage, referenceFolder, sampleId.ToString(), null,
                    samplesInfo[sampleId].Ploidy, QualityFilterThreshold, isPedigreeInfoSupplied, denovoQualityThreshold);
            }
            return 0;
        }

        private List<SampleList<CanvasSegment>> GetHighestLikelihoodSegments(List<SampleList<OverlappingSegmentsRegion>> segmentSetsFromCommonCnvs,
            SampleList<SamplesInfo> pedigreeMembersInfo, SampleList<CopyNumberModel> copyNumberModel)
        {

            Parallel.ForEach(
                segmentSetsFromCommonCnvs,
                new ParallelOptions
                {
                    MaxDegreeOfParallelism = Math.Min(Environment.ProcessorCount, CallerParameters.MaxCoreNumber)
                },
                segmentSet => GetHighestLikelihoodSegmentsSet(segmentSet, pedigreeMembersInfo, copyNumberModel)
            );
            return segmentSetsFromCommonCnvs.Select(sampleList => sampleList.Select(x => x.Value.GetSet().Select(y => (x.Key, y))).
                ZipMany(sampleRegion => sampleRegion.ToSampleList())).SelectMany(x => x).ToList();
        }
        /// <summary>
        /// Derives metrics from b-allele counts within each segment and determines whereas to use them for calculating MCC
        /// </summary>
        /// <param name="canvasSegments"></param>
        /// <param name="segmentIndex"></param>
        /// <returns></returns>
        private bool UseMafInformation(SampleList<CanvasSegment> canvasSegments)
        {
            var alleles = canvasSegments.SampleData.Select(segments => segments.Balleles?.TotalCoverage);
            var alleleCounts = alleles.Select(allele => allele?.Count ?? 0).ToList();
            bool lowAlleleCounts = alleleCounts.Select(x => x < CallerParameters.DefaultReadCountsThreshold).Any(c => c == true);
            var coverageCounts = canvasSegments.SampleData.Select(segments => segments.MedianCount).ToList();
            var isSkewedHetHomRatio = false;
            double alleleDensity = canvasSegments.SampleData.First().Length /
                                   Math.Max(alleleCounts.Average(), 1.0);
            bool useCnLikelihood = lowAlleleCounts ||
                                   alleleDensity < CallerParameters.DefaultAlleleDensityThreshold ||
                                   alleleCounts.Any(x => x > CallerParameters.DefaultPerSegmentAlleleMaxCounts) ||
                                   coverageCounts.Any(coverage => coverage < MedianCoverageThreshold) ||
                                   isSkewedHetHomRatio;
            // for now only use lowAlleleCounts metric
            return lowAlleleCounts;
        }


        private static SampleList<List<CanvasSegment>> MergeSegments(SampleList<List<CanvasSegment>> segments, int minimumCallSize)
        {
            int nSegments = segments.First().Value.Count;
            var copyNumbers = new List<List<int>>(nSegments);
            var qscores = new List<double>(nSegments);
            foreach (int segmentIndex in Enumerable.Range(0, nSegments))
            {
                copyNumbers.Add(segments.Select(s => s.Value[segmentIndex].CopyNumber).ToList());
                qscores.Add(segments.Select(s => s.Value[segmentIndex].QScore).Average());
            }

            if (copyNumbers == null && qscores != null || copyNumbers != null & qscores == null)
                throw new ArgumentException("Both copyNumbers and qscores arguments must be specified.");
            if (copyNumbers != null && copyNumbers.Count != nSegments)
                throw new ArgumentException("Length of copyNumbers list should be equal to the number of segments.");
            if (qscores != null && qscores.Count != nSegments)
                throw new ArgumentException("Length of qscores list should be equal to the number of segments.");

            var mergedSegments = new SampleList<List<CanvasSegment>>();
            foreach (var sampleSegments in segments)
            {
                var mergedAllSegments = CanvasSegment.MergeSegments(sampleSegments.Value.ToList(),
                    minimumCallSize, 10000, copyNumbers, qscores);
                mergedSegments.Add(sampleSegments.Key, mergedAllSegments);
            }
            return mergedSegments;
        }

        /// <summary>
        /// Create CanvasSegments from common CNVs bed file and overlap with CanvasPartition
        /// segments to create SegmentHaplotypes
        /// </summary>
        private List<SampleList<OverlappingSegmentsRegion>> CreateSegmentSetsFromCommonCnvs(SampleList<string> variantFrequencyFiles,
            int defaultAlleleCountThreshold, string commonCNVsbedPath, SampleList<Segments> sampleSegments)
        {
            if (commonCNVsbedPath == null)
            {
                var defaultSampleRegions = sampleSegments
                    .SelectData(segments => segments.AllSegments.Select(segment => new OverlappingSegmentsRegion(segment)).ToList());
                return GetOverlappingSegmentsRegionSampleLists(defaultSampleRegions);
            }

            var commonRegions = ReadCommonRegions(commonCNVsbedPath);
            var chromosomes = sampleSegments.SampleData.First().GetChromosomes();
            if (IsIdenticalChromosomeNames(commonRegions, chromosomes))
                throw new ArgumentException(
                    $"Chromosome names in a common CNVs bed file {commonCNVsbedPath} does not match the genome reference");

            var segmentIntervalsByChromosome = new Dictionary<string, List<BedInterval>>();
            var genomicBinsByChromosome = new Dictionary<string, IReadOnlyList<SampleGenomicBin>>();

            Parallel.ForEach(
                chromosomes,
                chr =>
                {
                    genomicBinsByChromosome[chr] = sampleSegments.SampleData.First().GetGenomicBinsForChromosome(chr);
                    segmentIntervalsByChromosome[chr] =
                        CanvasSegment.RemapGenomicToBinCoordinates(commonRegions[chr], genomicBinsByChromosome[chr]);
                });

            var sampleRegions = new SampleList<List<OverlappingSegmentsRegion>>();
            foreach (var sampleId in sampleSegments.SampleIds)
            {
                var commonIntervals = commonRegions.ToDictionary(kvp => kvp.Key, kvp => kvp.Value.Select(bedEntry => bedEntry.Interval).ToList());
                var allelesByChromosomeCommonSegs = CanvasIO.ReadFrequenciesWrapper(_logger,
                    new FileLocation(variantFrequencyFiles[sampleId]), commonIntervals);
                var segmentsSets = GetSegmentSets(defaultAlleleCountThreshold, commonRegions,
                    genomicBinsByChromosome, segmentIntervalsByChromosome, allelesByChromosomeCommonSegs, sampleSegments[sampleId]);
                sampleRegions.Add(sampleId, segmentsSets);
            }

            return GetOverlappingSegmentsRegionSampleLists(sampleRegions);
        }

        private static List<SampleList<OverlappingSegmentsRegion>> GetOverlappingSegmentsRegionSampleLists(SampleList<List<OverlappingSegmentsRegion>> sampleRegions)
        {
            return sampleRegions
                .SelectData((sampleId, regions) => regions.Select(region => (sampleId, region))).SampleData
                .ZipMany(sampleRegion => sampleRegion.ToSampleList())
                .ToList();
        }

        private static List<OverlappingSegmentsRegion> GetSegmentSets(int defaultAlleleCountThreshold, Dictionary<string, List<BedEntry>> commonRegions,
            Dictionary<string, IReadOnlyList<SampleGenomicBin>> genomicBinsByChromosome, Dictionary<string, List<BedInterval>> segmentIntervalsByChromosome,
            Dictionary<string, List<Balleles>> allelesByChromosomeCommonSegs, Segments segments)
        {
            var segmentsSetByChromosome = new ConcurrentDictionary<string, List<OverlappingSegmentsRegion>>();
            Parallel.ForEach(
                segments.GetChromosomes(),
                chr =>
                {
                    var segmentsByChromosome = segments.GetSegmentsForChromosome(chr).ToList();

                    if (commonRegions.Keys.Any(chromosome => chromosome == chr))
                    {
                        var commonCnvCanvasSegments = CanvasSegment.CreateSegmentsFromCommonCnvs(genomicBinsByChromosome[chr],
                            segmentIntervalsByChromosome[chr], allelesByChromosomeCommonSegs[chr]);

                        segmentsSetByChromosome[chr] = CanvasSegment.MergeCommonCnvSegments(segmentsByChromosome,
                            commonCnvCanvasSegments, defaultAlleleCountThreshold);
                    }
                    else
                    {
                        segmentsSetByChromosome[chr] = segmentsByChromosome.Select(
                            segment => new OverlappingSegmentsRegion(new List<CanvasSegment> { segment }, null)).ToList();
                    }
                });
            return segmentsSetByChromosome.OrderBy(i => i.Key).Select(x => x.Value).SelectMany(x => x).ToList();
        }

        private static Dictionary<string, List<BedEntry>> ReadCommonRegions(string commonCNVsbedPath)
        {
            Dictionary<string, List<BedEntry>> commonRegions;
            using (var reader = new BedReader(new GzipOrTextReader(commonCNVsbedPath)))
            {
                var commonTmpRegions = reader.LoadAllEntries();
                commonRegions = Utilities.SortAndOverlapCheck(commonTmpRegions, commonCNVsbedPath);
            }
            return commonRegions;
        }

        private static bool IsIdenticalChromosomeNames(Dictionary<string, List<BedEntry>> commonRegions,
            ICollection<string> chromsomeNames)
        {
            var chromsomes = new HashSet<string>(chromsomeNames);
            return commonRegions.Keys.Count(chromosome => chromsomes.Contains(chromosome)) == 0;
        }

        private void EstimateQScores(SampleList<CanvasSegment> canvasSegments, SampleList<SamplesInfo> pedigreeMembersInfo,
            PedigreeInfo pedigreeInfo, CopyNumbersLikelihoods copyNumberLikelihoods, SampleList<int> copyNumbers)
        {
            foreach (var sampleId in canvasSegments.SampleIds)
            {
                canvasSegments[sampleId].QScore = GetSingleSampleQualityScore(copyNumberLikelihoods.SingleSampleLikelihoods[sampleId], copyNumbers[sampleId], sampleId.ToString());
                canvasSegments[sampleId].CopyNumber = copyNumbers[sampleId];
                if (canvasSegments[sampleId].QScore < QualityFilterThreshold)
                    canvasSegments[sampleId].Filter = $"q{QualityFilterThreshold}";
            }
            if (pedigreeInfo != null)
                SetDenovoQualityScores(canvasSegments, pedigreeMembersInfo, pedigreeInfo.ParentsIds, pedigreeInfo.OffspringsIds, copyNumberLikelihoods);
        }

        private void SetDenovoQualityScores(SampleList<CanvasSegment> canvasSegments, SampleList<SamplesInfo> samplesInfo, List<SampleId> parentIDs, List<SampleId> offspringIDs,
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
                if (Double.IsInfinity(deNovoQualityScore) | deNovoQualityScore > CallerParameters.MaxQscore)
                    deNovoQualityScore = CallerParameters.MaxQscore;
                canvasSegments[probandId].DqScore = deNovoQualityScore;
            }
        }

        private bool IsPassVariant(SampleList<CanvasSegment> canvasSegments, SampleId sampleId)
        {
            return canvasSegments[sampleId].QScore > QualityFilterThreshold;
        }

        private bool IsCommonCnv(SampleList<CanvasSegment> canvasSegments, SampleList<SamplesInfo> samplesInfo, List<SampleId> parentIDs, SampleId probandId)
        {
            int parent1CopyNumber = GetCnState(canvasSegments, parentIDs.First(), CallerParameters.MaximumCopyNumber);
            int parent2CopyNumber = GetCnState(canvasSegments, parentIDs.Last(), CallerParameters.MaximumCopyNumber);
            int probandCopyNumber = GetCnState(canvasSegments, probandId, CallerParameters.MaximumCopyNumber);
            var parent1Genotypes = GenerateCnAlleles(parent1CopyNumber);
            var parent2Genotypes = GenerateCnAlleles(parent2CopyNumber);
            var probandGenotypes = GenerateCnAlleles(probandCopyNumber);
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

        private bool IsReferenceVariant(SampleList<CanvasSegment> canvasSegments, SampleList<SamplesInfo> samplesInfo, SampleId sampleId)
        {
            var segment = canvasSegments[sampleId];
            return GetCnState(canvasSegments, sampleId, CallerParameters.MaximumCopyNumber) == samplesInfo[sampleId].GetPloidy(segment);
        }


        private static int GetCnState(SampleList<CanvasSegment> canvasSegmentsSet, SampleId sampleId, int maximumCopyNumber)
        {
            return Math.Min(canvasSegmentsSet[sampleId].CopyNumber, maximumCopyNumber - 1);
        }

        public static int AggregateVariantCoverage(ref List<CanvasSegment> segments)
        {
            var variantCoverage = segments.SelectMany(segment => segment.Balleles.TotalCoverage).ToList();
            return variantCoverage.Any() ? Utilities.Median(variantCoverage) : 0;
        }


        /// <summary>
        /// Identify variant with the highest likelihood at a given setPosition and assign relevant scores
        /// </summary>
        private void GetHighestLikelihoodSegmentsSet(SampleList<OverlappingSegmentsRegion> canvasSegmentsSet, SampleList<SamplesInfo> pedigreeMembersInfo,
            SampleList<CopyNumberModel> model)
        {
            SegmentsSet segmentSet;

            if (canvasSegmentsSet.SampleData.First().SetA == null)
                segmentSet = SegmentsSet.SetB;
            else if (canvasSegmentsSet.SampleData.First().SetB == null)
                segmentSet = SegmentsSet.SetA;
            else
                segmentSet = GetSegmentSetLikelihood(canvasSegmentsSet, pedigreeMembersInfo, model,
                                 SegmentsSet.SetA) >
                             GetSegmentSetLikelihood(canvasSegmentsSet, pedigreeMembersInfo, model,
                                 SegmentsSet.SetB)
                    ? SegmentsSet.SetA
                    : SegmentsSet.SetB;

            canvasSegmentsSet.SampleIds.ForEach(id => canvasSegmentsSet[id].SetSet(segmentSet));
        }

        private double GetSegmentSetLikelihood(SampleList<OverlappingSegmentsRegion> canvasSegmentsSet, SampleList<SamplesInfo> samplesInfo,
            SampleList<CopyNumberModel> copyNumberModel, SegmentsSet segmentsSet)
        {
            double segmentSetLikelihood = 0;
            foreach (var sampleId in canvasSegmentsSet.SampleIds)
                canvasSegmentsSet[sampleId].SetSet(segmentsSet);

            var canvasSegments = new List<SampleList<CanvasSegment>>();
            int nSegments = canvasSegmentsSet.First().Value.GetSet().Count;
            for (var canvasSegmentIndex = 0; canvasSegmentIndex < nSegments; canvasSegmentIndex++)
            {
                var canvasSegment = new SampleList<CanvasSegment>();
                foreach (var id in canvasSegmentsSet.SampleIds)
                    canvasSegment.Add(id, canvasSegmentsSet[id].GetSet()[canvasSegmentIndex]);
                canvasSegments.Add(canvasSegment);
            }
            foreach (var canvasSegment in canvasSegments)
            {
                var copyNumbersLikelihoods = GetCopyNumbersLikelihoods(canvasSegment, samplesInfo, copyNumberModel);
                var copyNumbers = GetCopyNumbersNoPedigreeInfo(canvasSegment, copyNumbersLikelihoods);
                segmentSetLikelihood += copyNumbersLikelihoods.MaximalLikelihood;
            }

            segmentSetLikelihood /= nSegments;
            return segmentSetLikelihood;
        }

        /// <summary>
        /// Identify variant with the highest likelihood at a given setPosition and assign relevant scores
        /// </summary>
        public void CallVariant(SampleList<CanvasSegment> canvasSegments, SampleList<SamplesInfo> samplesInfo,
            SampleList<CopyNumberModel> copyNumberModel,
            PedigreeInfo pedigreeInfo, Dictionary<int, List<Genotype>> genotypes)
        {
            var copyNumbersLikelihoods = GetCopyNumbersLikelihoods(canvasSegments, samplesInfo, copyNumberModel);

            var copyNumbers = pedigreeInfo != null
                ? GetCopyNumbersWithPedigreeInfo(canvasSegments, copyNumberModel, pedigreeInfo, copyNumbersLikelihoods)
                : GetCopyNumbersNoPedigreeInfo(canvasSegments, copyNumbersLikelihoods);

            EstimateQScores(canvasSegments, samplesInfo, pedigreeInfo, copyNumbersLikelihoods, copyNumbers);

            // TODO: this will be integrated with GetCopyNumbers* on a model level as a part of https://jira.illumina.com/browse/CANV-404
            if (!UseMafInformation(canvasSegments) && pedigreeInfo != null)
                AssignMccWithPedigreeInfo(canvasSegments, samplesInfo, copyNumberModel, pedigreeInfo, genotypes);
            if (!UseMafInformation(canvasSegments) && pedigreeInfo == null)
                AssignMccNoPedigreeInfo(canvasSegments, samplesInfo, copyNumberModel, genotypes);                         

        }


        /// <summary>
        /// Calculates maximal likelihood for copy numbers. Updated CanvasSegment CopyNumber only. 
        /// </summary>
        public SampleList<int> GetCopyNumbersWithPedigreeInfo(SampleList<CanvasSegment> segments, SampleList<CopyNumberModel> model,
            PedigreeInfo pedigreeInfo, CopyNumbersLikelihoods copyNumbersLikelihoods)
        {
            var sampleCopyNumbers = new SampleList<int>();
            segments.SampleIds.ForEach(id => sampleCopyNumbers.Add(id, 2));
            var parent1Likelihood = copyNumbersLikelihoods.SingleSampleLikelihoods[pedigreeInfo.ParentsIds.First()];
            var parent2Likelihood = copyNumbersLikelihoods.SingleSampleLikelihoods[pedigreeInfo.ParentsIds.Last()];
            var copyNumbersRange = Enumerable.Range(0, CallerParameters.MaximumCopyNumber).ToList();
            foreach (int copyNumberParent1 in copyNumbersRange)
            {
                foreach (int copyNumberParent2 in copyNumbersRange)
                {
                    foreach (var offspringGtStates in pedigreeInfo.OffspringsGenotypes)
                    {
                        double currentLikelihood = parent1Likelihood[copyNumberParent1] * parent2Likelihood[copyNumberParent2];
                        for (var counter = 0; counter < pedigreeInfo.OffspringsIds.Count; counter++)
                        {
                            var child = pedigreeInfo.OffspringsIds[counter];
                            int copyNumberChild = Math.Min(offspringGtStates[counter].CountsA + offspringGtStates[counter].CountsB,
                                    CallerParameters.MaximumCopyNumber - 1);
                            currentLikelihood *= pedigreeInfo.TransitionMatrix[copyNumberParent1][offspringGtStates[counter].CountsA] *
                                                 pedigreeInfo.TransitionMatrix[copyNumberParent2][offspringGtStates[counter].CountsB] *
                                                 copyNumbersLikelihoods.SingleSampleLikelihoods[child][copyNumberChild];
                        }

                        currentLikelihood = Double.IsNaN(currentLikelihood) || Double.IsInfinity(currentLikelihood)
                            ? 0
                            : currentLikelihood;
                        var copyNumbers = new[] { copyNumberParent1, copyNumberParent2 }
                            .Concat(offspringGtStates.Select(x => x.CountsA + x.CountsB)).ToArray();
                        copyNumbersLikelihoods.SetJointProbability(currentLikelihood, copyNumbers);

                        if (currentLikelihood > copyNumbersLikelihoods.MaximalLikelihood)
                        {
                            sampleCopyNumbers = new SampleList<int>();
                            copyNumbersLikelihoods.MaximalLikelihood = currentLikelihood;
                            sampleCopyNumbers.Add(pedigreeInfo.ParentsIds.First(), copyNumberParent1);
                            sampleCopyNumbers.Add(pedigreeInfo.ParentsIds.Last(), copyNumberParent2);
                            for (int counter = 0; counter < pedigreeInfo.OffspringsIds.Count; counter++)
                            {
                                sampleCopyNumbers.Add(pedigreeInfo.OffspringsIds[counter], offspringGtStates[counter].CountsA + offspringGtStates[counter].CountsB);
                            }
                        }
                    }
                }
            }
            return sampleCopyNumbers;
        }


        /// <summary>
        /// Calculates maximal likelihood for copy numbers. Updated CanvasSegment CopyNumber only. 
        /// </summary>
        public SampleList<int> GetCopyNumbersNoPedigreeInfo(SampleList<CanvasSegment> segments, CopyNumbersLikelihoods copyNumbersLikelihoods)
        {
            var sampleCopyNumbers = new SampleList<int>();
            double maximalLikelihood = 1;
            foreach (var sampleId in segments.SampleIds)
            {
                double maxSampleLikelihoods = copyNumbersLikelihoods.SingleSampleLikelihoods[sampleId].Max(x => x.Value);
                maximalLikelihood *= maxSampleLikelihoods;
                int copyNumber = copyNumbersLikelihoods.SingleSampleLikelihoods[sampleId].First(x => x.Value == maxSampleLikelihoods).Key;
                sampleCopyNumbers.Add(sampleId, copyNumber);
            }
            copyNumbersLikelihoods.MaximalLikelihood = maximalLikelihood;
            return sampleCopyNumbers;
        }


        /// <summary>
        /// Calculates maximal likelihood for genotypes given a copy number call. Updated MajorChromosomeCount.
        /// </summary>
        public void AssignMccWithPedigreeInfo(SampleList<CanvasSegment> canvasSegments,
            SampleList<SamplesInfo> samplesInfo, SampleList<CopyNumberModel> model, PedigreeInfo pedigreeInfo,
            Dictionary<int, List<Genotype>> genotypes)
        {
            double maximalLikelihood = Double.MinValue;
            int parent1CopyNumber = canvasSegments[pedigreeInfo.ParentsIds.First()].CopyNumber;
            int parent2CopyNumber = canvasSegments[pedigreeInfo.ParentsIds.Last()].CopyNumber;

            foreach (var parent1GtStates in genotypes[parent1CopyNumber])
            {
                foreach (var parent2GtStates in genotypes[parent2CopyNumber])
                {
                    var bestChildGtStates = new List<Genotype>();
                    double currentLikelihood = 1;
                    foreach (SampleId child in pedigreeInfo.OffspringsIds)
                    {
                        int childCopyNumber = canvasSegments[child].CopyNumber;
                        bool isInheritedCnv = !canvasSegments[child].DqScore.HasValue;
                        double bestLikelihood = Double.MinValue;
                        Genotype bestGtState = null;
                        bestLikelihood = GetProbandLikelihood(model[child], genotypes, childCopyNumber,
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
                        AssignMCC(canvasSegments[pedigreeInfo.ParentsIds.First()], model[pedigreeInfo.ParentsIds.First()], genotypes, parent1GtStates, parent1CopyNumber);
                        AssignMCC(canvasSegments[pedigreeInfo.ParentsIds.Last()], model[pedigreeInfo.ParentsIds.Last()], genotypes, parent2GtStates, parent2CopyNumber);
                        var counter = 0;
                        foreach (SampleId child in pedigreeInfo.OffspringsIds)
                        {
                            if (bestChildGtStates[counter] == null) continue;
                            int childCopyNumber = canvasSegments[child].CopyNumber;
                            AssignMCC(canvasSegments[child], model[child], genotypes, bestChildGtStates[counter], childCopyNumber);
                            counter++;
                        }
                    }
                }
            }
        }

        private double GetProbandLikelihood(CopyNumberModel copyNumberModel, Dictionary<int, List<Genotype>> genotypes,
            int childCopyNumber, Genotype parent1GtStates, Genotype parent2GtStates, bool isInheritedCnv, CanvasSegment canvasSegment,
            double bestLikelihood, ref Genotype bestGtState)
        {
            foreach (var childGtState in genotypes[childCopyNumber])
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

        private static void AssignMCC(CanvasSegment canvasSegment, CopyNumberModel copyNumberModel,
            Dictionary<int, List<Genotype>> genotypes, Genotype gtStates, int copyNumber)
        {
            const int diploidCopyNumber = 2;
            const int haploidCopyNumber = 1;
            if (copyNumber > diploidCopyNumber)
            {

                canvasSegment.MajorChromosomeCount =
                    Math.Max(gtStates.CountsA, gtStates.CountsB);
                int? selectedGtState = genotypes[copyNumber].IndexOf(gtStates);
                canvasSegment.MajorChromosomeCountScore =
                    copyNumberModel.GetGtLikelihoodScore(canvasSegment.Balleles.GetAlleleCounts(), genotypes[copyNumber], ref selectedGtState);
                copyNumberModel.GetGtLikelihoodScore(canvasSegment.Balleles.GetAlleleCounts(), genotypes[copyNumber], ref selectedGtState);
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

        public bool IsGtPedigreeConsistent(Genotype parentGtStates, Genotype childGtStates)
        {
            if (parentGtStates.CountsA == childGtStates.CountsA || parentGtStates.CountsB == childGtStates.CountsA ||
                parentGtStates.CountsA == childGtStates.CountsB || parentGtStates.CountsB == childGtStates.CountsB)
                return true;
            return false;
        }

        /// <summary>
        /// Calculates maximal likelihood for segments without SNV allele ratios. Updated CanvasSegment CopyNumber only. 
        /// </summary>
        public CopyNumbersLikelihoods GetCopyNumbersLikelihoods(SampleList<CanvasSegment> canvasSegments, SampleList<SamplesInfo> samplesInfo,
            SampleList<CopyNumberModel> copyNumberModel)
        {
            const double maxCoverageMultiplier = 3.0;
            var singleSampleLikelihoods = new SampleList<Dictionary<int, double>>();
            foreach (var sampleId in canvasSegments.SampleIds)
            {
                var density = new Dictionary<int, double>();
                foreach (int copyNumber in Enumerable.Range(0, CallerParameters.MaximumCopyNumber))
                {
                    double currentLikelihood =
                        copyNumberModel[sampleId].GetCnLikelihood(
                            Math.Min(canvasSegments[sampleId].MedianCount,
                                samplesInfo[sampleId].MeanCoverage * maxCoverageMultiplier))[copyNumber];
                    currentLikelihood = Double.IsNaN(currentLikelihood) || Double.IsInfinity(currentLikelihood)
                        ? 0
                        : currentLikelihood;
                    density[copyNumber] = currentLikelihood;
                }
                singleSampleLikelihoods.Add(sampleId, density);
            }
            return new CopyNumbersLikelihoods(singleSampleLikelihoods, CallerParameters.MaximumCopyNumber);
        }

        /// <summary>
        /// Calculates maximal likelihood for segments with SNV allele counts given CopyNumber. Updated MajorChromosomeCount.
        /// </summary>   
        public void AssignMccNoPedigreeInfo(SampleList<CanvasSegment> canvasSegments, SampleList<SamplesInfo> pedigreeMembersInfo,
            SampleList<CopyNumberModel> model, Dictionary<int, List<Genotype>> genotypes)
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
        /// Generate all possible copy number combinations with the maximal number of copy numbers per segment set to maxAlleleNumber.
        /// </summary>
        /// <param name="numberOfCnStates"></param>
        /// <param name="maxAlleleNumber"></param>
        /// <returns></returns>
        public static List<List<int>> GenerateCopyNumberCombinations(int numberOfCnStates, int maxAlleleNumber)
        {
            if (numberOfCnStates <= 0)
                throw new ArgumentOutOfRangeException(nameof(numberOfCnStates));
            var cnStates = Enumerable.Range(0, numberOfCnStates).ToList();
            var allCombinations = new List<List<int>>();
            for (int currentAlleleNumber = 1; currentAlleleNumber <= maxAlleleNumber; currentAlleleNumber++)
            {
                var currentCombination = new Combinations<int>(cnStates, currentAlleleNumber);
                var list = currentCombination.Select(x => x.ToList()).ToList();
                allCombinations.AddRange(list);
            }
            return allCombinations;
        }


        /// <summary>
        /// Generate all possible copy number genotype combinations with the maximal number of alleles per segment set to maxAlleleNumber.
        /// </summary>
        /// <param name="copyNumber"></param>
        /// <returns> </returns>
        public static List<int> GenerateCnAlleles(int copyNumber)
        {
            if (copyNumber == 0)
                return new List<int> { 0 };

            if (copyNumber == 1)
                return new List<int> { 0, 1 };

            var alleles = new List<int>();
            for (int allele = 1; allele <= copyNumber; allele++)
                alleles.Add(allele);

            return alleles;
        }

        /// <summary>
        /// Generate all possible copy number genotype combinations with the maximal number of alleles per segment set to maxAlleleNumber.
        /// </summary>
        /// <param name="numberOfCnStates"></param>
        /// <returns> </returns>
        public Dictionary<int, List<Genotype>> GenerateGenotypeCombinations(int numberOfCnStates)
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


        public double GetTransitionProbability(int gt1Parent, int gt2Parent, int gt1Offspring, int gt2Offspring)
        {
            if (gt1Parent == gt1Offspring || gt1Parent == gt2Offspring ||
                gt2Parent == gt1Offspring || gt2Parent == gt2Offspring)
                return 0.5;
            return CallerParameters.DeNovoRate;
        }

        public enum Kinship
        {
            Other = 0,
            Parent = 1,
            Proband = 2
        }
        public static SampleList<Kinship> ReadPedigreeFile(string pedigreeFile)
        {
            var kinships = new SampleList<Kinship>();
            using (FileStream stream = new FileStream(pedigreeFile, FileMode.Open, FileAccess.Read))
            using (StreamReader reader = new StreamReader(stream))
            {
                string row;
                while ((row = reader.ReadLine()) != null)
                {
                    string[] fields = row.Split('\t');
                    string maternalId = fields[2];
                    string paternallId = fields[3];
                    string proband = fields[5];
                    if (maternalId == "0" && paternallId == "0")
                        kinships.Add(new SampleId(fields[1]), Kinship.Parent);
                    else if (proband == "affected")
                        kinships.Add(new SampleId(fields[1]), Kinship.Proband);
                    else
                        Console.WriteLine($"Unused pedigree member: {row}");
                }
            }
            return kinships;
        }

        public double GetSingleSampleQualityScore(Dictionary<int, double> copyNumbersLikelihoods, int cnState, string sampleName)
        {
            double normalizationConstant = copyNumbersLikelihoods.Select(ll => ll.Value).Sum();
            double qscore = -10.0 * Math.Log10((normalizationConstant - copyNumbersLikelihoods[cnState]) / normalizationConstant);
            if (Double.IsInfinity(qscore) | qscore > CallerParameters.MaxQscore)
                qscore = CallerParameters.MaxQscore;
            return qscore;
        }

        public double GetConditionalDeNovoQualityScore(CopyNumbersLikelihoods copyNumbersLikelihoods, SampleId probandId, SampleList<CanvasSegment> canvasSegments,
            List<SampleId> parentIDs)
        {

            var numerator = 0.0;
            var denominator = 0.0;
            const int diploidState = 2;
            var names = copyNumbersLikelihoods.SampleNames;
            int parent1Index = names.IndexOf(parentIDs.First().ToString());
            int parent2Index = names.IndexOf(parentIDs.Last().ToString());
            int probandIndex = names.IndexOf(probandId.ToString());

            var probandMarginalProbabilities = copyNumbersLikelihoods.GetMarginalProbability(CallerParameters.MaximumCopyNumber, probandId.ToString());
            int probandCopyNumber = GetCnState(canvasSegments, probandId, CallerParameters.MaximumCopyNumber);
            double normalization = probandMarginalProbabilities[probandCopyNumber] + probandMarginalProbabilities[diploidState];
            double probandMarginalAlt = probandMarginalProbabilities[probandCopyNumber] / normalization;

            foreach (var copyNumberIndex in copyNumbersLikelihoods.Indices.Where(x => x[probandIndex] == probandCopyNumber))
            {
                if (!(copyNumbersLikelihoods.GetJointProbability(copyNumberIndex.ToArray()) > 0.0))
                    continue;
                double holder = copyNumbersLikelihoods.GetJointProbability(copyNumberIndex);
                denominator += holder;
                if (copyNumberIndex[parent1Index] == diploidState && copyNumberIndex[parent2Index] == diploidState)
                    numerator += holder;
            }

            const double q60 = 0.000001;
            double denovoProbability = (1 - numerator / denominator) * (1 - probandMarginalAlt);
            return -10.0 * Math.Log10(Math.Max(denovoProbability, q60));
        }
    }
}



