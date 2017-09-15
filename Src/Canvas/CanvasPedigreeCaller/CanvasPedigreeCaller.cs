using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MathNet.Numerics.Distributions;
using System.Threading.Tasks;
using CanvasCommon;
using Combinatorics.Collections;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Isas.Framework.DataTypes;
using Isas.Framework.Logging;
using Isas.SequencingFiles;
using Isas.SequencingFiles.Vcf;
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

        internal int CallVariantsInPedigree(List<string> variantFrequencyFiles, List<string> segmentFiles,
            string outVcfFile, string ploidyBedPath, string referenceFolder, List<string> sampleNames, string commonCNVsbedPath, string pedigreeFile)
        {
            // load files
            // initialize data structures and classes
            int fileCounter = 0;
            var kinships = ReadPedigreeFile(pedigreeFile);
            var pedigreeMembersInfo = new SampleList<PedigreeMemberInfo>();
            var sampleSegments = new SampleList<Segments>();
            var copyNumberModels = new SampleList<CopyNumberModel>();
            var variantFrequencyFilesSampleList = new SampleList<string>();


            foreach (string sampleName in sampleNames)
            {
                var sampleId = new SampleId(sampleName);
                var segment = Segments.ReadSegments(_logger, new FileLocation(segmentFiles[fileCounter]));
                segment.AddAlleles(CanvasIO.ReadFrequenciesWrapper(_logger, new FileLocation(variantFrequencyFiles[fileCounter]), segment.IntervalsByChromosome));
                sampleSegments.Add(sampleId, segment);
                var pedigreeMemberInfo = PedigreeMemberInfo.GetPedigreeMemberInfo(segment, ploidyBedPath, CallerParameters.NumberOfTrimmedBins, sampleId);
                var copyNumberModel = new CopyNumberModel(CallerParameters.MaximumCopyNumber, pedigreeMemberInfo);
                pedigreeMembersInfo.Add(sampleId, pedigreeMemberInfo);
                copyNumberModels.Add(sampleId, copyNumberModel);
                variantFrequencyFilesSampleList.Add(sampleId, variantFrequencyFiles[fileCounter]);
                fileCounter++;
            }
            var segmentSetsFromCommonCnvs = CreateSegmentSetsFromCommonCnvs(variantFrequencyFilesSampleList,
                CallerParameters.DefaultReadCountsThreshold, commonCNVsbedPath, sampleSegments);

            var parentsIds = kinships.Where(kin => kin.Value.Equals(Kinship.Parent)).Select(kin => kin.Key).ToList();
            var offspringsIds = kinships.Where(kin => kin.Value.Equals(Kinship.Proband)).Select(kin => kin.Key).ToList();
            var parentalGenotypes = GenerateParentalGenotypes(CallerParameters.MaximumCopyNumber);
            var offspringsGenotypes = new List<List<Genotype>>(Convert.ToInt32(Math.Pow(parentalGenotypes.Count, offspringsIds.Count)));
            GenerateOffspringGenotypes(offspringsGenotypes, parentalGenotypes, offspringsIds.Count, new List<Genotype>());
            var genotypes = GenerateGenotypeCombinations(CallerParameters.MaximumCopyNumber);

            if (offspringsGenotypes.Count > CallerParameters.MaxNumOffspringGenotypes)
            {
                offspringsGenotypes.Shuffle();
                offspringsGenotypes = offspringsGenotypes.Take(CallerParameters.MaxNumOffspringGenotypes).ToList();
            }

            double[][] transitionMatrix = GetTransitionMatrix(CallerParameters.MaximumCopyNumber);
            Parallel.ForEach(
                segmentSetsFromCommonCnvs,
                new ParallelOptions
                {
                    MaxDegreeOfParallelism = Math.Min(Environment.ProcessorCount, CallerParameters.MaxCoreNumber)
                },
                segmentSet => CallVariantInPedigree(segmentSet, pedigreeMembersInfo, copyNumberModels, parentsIds, offspringsIds,
                    transitionMatrix, offspringsGenotypes, genotypes)
            );

            var variantCalledSegments = new SampleList<List<CanvasSegment>>();
            foreach (var key in pedigreeMembersInfo.SampleIds)
                variantCalledSegments.Add(key, segmentSetsFromCommonCnvs.SelectMany(segmentsSet => segmentsSet[key].GetSet()).ToList());

            var mergedVariantCalledSegments = MergeSegments(variantCalledSegments, CallerParameters.MinimumCallSize);
            var outputFolder = new FileLocation(outVcfFile).Directory;
            foreach (var sampleId in offspringsIds.Union(parentsIds))
            {
                var coverageOutputPath = SingleSampleCallset.GetCoverageAndVariantFrequencyOutput(outputFolder,
                    sampleId.ToString());
                CanvasSegment.WriteCoveragePlotData(mergedVariantCalledSegments[sampleId], pedigreeMembersInfo[sampleId].MeanCoverage,
                    pedigreeMembersInfo[sampleId].Ploidy, coverageOutputPath, referenceFolder);
            }

            var ploidies = offspringsIds.Union(parentsIds).Select(id => pedigreeMembersInfo[id].Ploidy).ToList();
            var diploidCoverage = offspringsIds.Union(parentsIds).Select(id => pedigreeMembersInfo[id].MeanCoverage).ToList();
            var names = offspringsIds.Union(parentsIds).Select(x => x.ToString()).ToList();
            CanvasSegmentWriter.WriteMultiSampleSegments(outVcfFile, mergedVariantCalledSegments, diploidCoverage, referenceFolder, names,
                null, ploidies, QualityFilterThreshold, isPedigreeInfoSupplied: true, denovoQualityThreshold: DeNovoQualityFilterThreshold);

            outputFolder = new FileLocation(outVcfFile).Directory;
            foreach (var sampleId in offspringsIds.Union(parentsIds))
            {
                var outputVcfPath = SingleSampleCallset.GetSingleSamplePedigreeVcfOutput(outputFolder, sampleId.ToString());
                CanvasSegmentWriter.WriteSegments(outputVcfPath.FullName, mergedVariantCalledSegments[sampleId],
                    pedigreeMembersInfo[sampleId].MeanCoverage, referenceFolder, sampleId.ToString(), null,
                    pedigreeMembersInfo[sampleId].Ploidy, QualityFilterThreshold,
                    isPedigreeInfoSupplied: true, denovoQualityThreshold: DeNovoQualityFilterThreshold);
            }
            return 0;
        }

        /// <summary>
        /// Derives metrics from b-allele counts within each segment and determines whereas to use them for calculating MCC
        /// </summary>
        /// <param name="canvasSegmentsSet"></param>
        /// <param name="segmentIndex"></param>
        /// <returns></returns>
        private bool UseMafInformation(SampleList<CanvasSegmentsSet> canvasSegmentsSet, int segmentIndex)
        {
            var alleles = canvasSegmentsSet.SampleData.Select(set => set.GetSet()[segmentIndex].Balleles?.TotalCoverage);
            var alleleCounts = alleles.Select(allele => allele?.Count ?? 0).ToList();
            bool lowAlleleCounts = alleleCounts.Select(x => x < CallerParameters.DefaultReadCountsThreshold).Any(c => c == true);
            var coverageCounts = canvasSegmentsSet.SampleData.Select(set => set.GetSet()[segmentIndex].MedianCount).ToList();
            var isSkewedHetHomRatio = false;
            double alleleDensity = canvasSegmentsSet.SampleData.First().GetSet()[segmentIndex].Length /
                                   Math.Max(alleleCounts.Average(), 1.0);
            bool useCnLikelihood = lowAlleleCounts ||
                                   alleleDensity < CallerParameters.DefaultAlleleDensityThreshold ||
                                   alleleCounts.Any(x => x > CallerParameters.DefaultPerSegmentAlleleMaxCounts) ||
                                   coverageCounts.Any(coverage => coverage < MedianCoverageThreshold) ||
                                   isSkewedHetHomRatio;
            // for now only use lowAlleleCounts metric
            return lowAlleleCounts;
        }


        internal int CallVariants(List<string> variantFrequencyFiles, List<string> segmentFiles, string outVcfFile,
            string ploidyBedPath, string referenceFolder, List<string> sampleNames, string commonCNVsbedPath)
        {
            // load files
            // initialize data structures and classes
            var fileCounter = 0;
            var pedigreeMembersInfo = new SampleList<PedigreeMemberInfo>();
            var sampleSegments = new SampleList<Segments>();
            var copyNumberModels = new SampleList<CopyNumberModel>();
            var variantFrequencyFilesSampleList = new SampleList<string>();


            foreach (string sampleName in sampleNames)
            {
                var sampleId = new SampleId(sampleName);
                var segment = Segments.ReadSegments(_logger, new FileLocation(segmentFiles[fileCounter]));
                segment.AddAlleles(CanvasIO.ReadFrequenciesWrapper(_logger, new FileLocation(variantFrequencyFiles[fileCounter]), segment.IntervalsByChromosome));
                sampleSegments.Add(sampleId, segment);
                var pedigreeMemberInfo = PedigreeMemberInfo.GetPedigreeMemberInfo(segment, ploidyBedPath, CallerParameters.NumberOfTrimmedBins, sampleId);
                var copyNumberModel = new CopyNumberModel(CallerParameters.MaximumCopyNumber, pedigreeMemberInfo);
                pedigreeMembersInfo.Add(sampleId, pedigreeMemberInfo);
                copyNumberModels.Add(sampleId, copyNumberModel);
                variantFrequencyFilesSampleList.Add(sampleId, variantFrequencyFiles[fileCounter]);
                fileCounter++;
            }
            var segmentSetsFromCommonCnvs = CreateSegmentSetsFromCommonCnvs(variantFrequencyFilesSampleList,
                CallerParameters.DefaultReadCountsThreshold, commonCNVsbedPath, sampleSegments);

            var genotypes = GenerateGenotypeCombinations(CallerParameters.MaximumCopyNumber);
            int maxAlleleNumber = Math.Min(CallerParameters.MaxAlleleNumber, pedigreeMembersInfo.Count());
            var copyNumberCombinations = GenerateCopyNumberCombinations(CallerParameters.MaximumCopyNumber, maxAlleleNumber);

            Parallel.ForEach(
                segmentSetsFromCommonCnvs,
                new ParallelOptions
                {
                    MaxDegreeOfParallelism = Math.Min(Environment.ProcessorCount, CallerParameters.MaxCoreNumber)
                },
                segmentSet => CallVariant(segmentSet, pedigreeMembersInfo, copyNumberModels, copyNumberCombinations, genotypes)
             );


            var variantCalledSegments = new SampleList<List<CanvasSegment>>();
            foreach (var key in pedigreeMembersInfo.SampleIds)
                variantCalledSegments.Add(key, segmentSetsFromCommonCnvs.SelectMany(segmentsSet => segmentsSet[key].GetSet()).ToList());

            var mergedVariantCalledSegments = MergeSegments(variantCalledSegments, CallerParameters.MinimumCallSize);
            var outputFolder = new FileLocation(outVcfFile).Directory;
            foreach (string sampleName in sampleNames)
            {
                var sampleId = new SampleId(sampleName);
                var coverageOutputPath = SingleSampleCallset.GetCoverageAndVariantFrequencyOutput(outputFolder, sampleName);
                CanvasSegment.WriteCoveragePlotData(mergedVariantCalledSegments[sampleId], pedigreeMembersInfo[sampleId].MeanCoverage,
                    pedigreeMembersInfo[sampleId].Ploidy, coverageOutputPath, referenceFolder);
                var outputVcfPath = SingleSampleCallset.GetSingleSamplePedigreeVcfOutput(outputFolder, sampleName);
                CanvasSegmentWriter.WriteSegments(outputVcfPath.FullName, mergedVariantCalledSegments[sampleId],
                    pedigreeMembersInfo[sampleId].MeanCoverage, referenceFolder, sampleName, null, pedigreeMembersInfo[sampleId].Ploidy,
                    QualityFilterThreshold, isPedigreeInfoSupplied: false);
            }
            return 0;
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
        /// <param name="variantFrequencyFiles"></param>
        /// <param name="segmentFile"></param>
        /// <param name="defaultAlleleCountThreshold"></param>
        /// <param name="commonCNVsbedPath"></param>
        /// <param name="pedigreeMember"></param>
        /// <returns></returns>
        private List<SampleList<CanvasSegmentsSet>> CreateSegmentSetsFromCommonCnvs(SampleList<string> variantFrequencyFiles,
            int defaultAlleleCountThreshold, string commonCNVsbedPath, SampleList<Segments> segments)
        {
            var canvasSegmentsSetBySample = new SampleList<List<CanvasSegmentsSet>>();

            var commonRegions = ReadCommonRegions(commonCNVsbedPath);
            var chromosomes = segments.SampleData.First().GetChromosomes();
            if (IsIdenticalChromosomeNames(commonRegions, chromosomes))
                throw new ArgumentException(
                    $"Chromosome names in a common CNVs bed file {commonCNVsbedPath} does not match the genome reference");

            var segmentIntervalsByChromosome = new Dictionary<string, List<BedInterval>>();
            var genomicBinsByChromosome = new Dictionary<string, IReadOnlyList<SampleGenomicBin>>();

            Parallel.ForEach(
                chromosomes,
                chr =>
                {
                    genomicBinsByChromosome[chr] = segments.SampleData.First().GetGenomicBinsForChromosome(chr);
                    segmentIntervalsByChromosome[chr] =
                        CanvasSegment.RemapGenomicToBinCoordinates(commonRegions[chr], genomicBinsByChromosome[chr]);
                });

            foreach (var sampleId in segments.SampleIds)
            {
                var allelesByChromosomeCommonSegs = CanvasIO.ReadFrequenciesWrapper(_logger,
                    new FileLocation(variantFrequencyFiles[sampleId]), segmentIntervalsByChromosome);
                var segmentsSetByChromosome = GetSegmentSets(defaultAlleleCountThreshold, commonRegions,
                    genomicBinsByChromosome, segmentIntervalsByChromosome, allelesByChromosomeCommonSegs, segments[sampleId]);
                canvasSegmentsSetBySample.Add(sampleId, segmentsSetByChromosome);
            }

            return GetSegmentsSetBySeampleId(canvasSegmentsSetBySample);
        }

        private static List<SampleList<CanvasSegmentsSet>> GetSegmentsSetBySeampleId(SampleList<List<CanvasSegmentsSet>> canvasSegmentsSetBySample)
        {
            int size = canvasSegmentsSetBySample.First().Value.Count;
            var canvasSegmentsSetBySegment = new List<SampleList<CanvasSegmentsSet>>(size);
            for (int segmentIndex = 0; segmentIndex < size; segmentIndex++)
            {
                canvasSegmentsSetBySegment.Add(new SampleList<CanvasSegmentsSet>());
            }
            foreach (var sampleId in canvasSegmentsSetBySample.SampleIds)
            {
                for (int segmentIndex = 0; segmentIndex < size; segmentIndex++)
                {
                    canvasSegmentsSetBySegment[segmentIndex].Add(sampleId, canvasSegmentsSetBySample[sampleId][segmentIndex]);
                }
            }
            return canvasSegmentsSetBySegment;
        }

        private static List<CanvasSegmentsSet> GetSegmentSets(int defaultAlleleCountThreshold, Dictionary<string, List<BedEntry>> commonRegions,
            Dictionary<string, IReadOnlyList<SampleGenomicBin>> genomicBinsByChromosome, Dictionary<string, List<BedInterval>> segmentIntervalsByChromosome,
            Dictionary<string, List<Balleles>> allelesByChromosomeCommonSegs, Segments segments)
        {
            var segmentsSetByChromosome = new ConcurrentDictionary<string, List<CanvasSegmentsSet>>();
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
                            segment => new CanvasSegmentsSet(new List<CanvasSegment> { segment }, null)).ToList();
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

        private void EstimateQScoresWithPedigreeInfo(SampleList<CanvasSegmentsSet> canvasSegmentsSet, SampleList<PedigreeMemberInfo> pedigreeMembersInfo,
            List<SampleId> parentIDs, List<SampleId> offspringIDs, int canvasSegmentIndex, CopyNumbersLikelihood copyNumberLikelihoods)
        {
            var sampleIds = parentIDs.Union(offspringIDs).ToList();
            foreach (var sampleId in sampleIds)
            {
                int cnState = GetCnState(canvasSegmentsSet, sampleId, canvasSegmentIndex, CallerParameters.MaximumCopyNumber);
                double singleSampleQualityScore = GetSingleSampleQualityScore(copyNumberLikelihoods, cnState, sampleId.ToString());
                canvasSegmentsSet[sampleId].GetSet()[canvasSegmentIndex].QScore = singleSampleQualityScore;
                if (canvasSegmentsSet[sampleId].GetSet()[canvasSegmentIndex].QScore < QualityFilterThreshold)
                    canvasSegmentsSet[sampleId].GetSet()[canvasSegmentIndex].Filter = $"q{QualityFilterThreshold}";
            }

            SetDenovoQualityScores(canvasSegmentsSet, pedigreeMembersInfo, parentIDs, offspringIDs, canvasSegmentIndex,
                copyNumberLikelihoods, sampleIds.Select(x => x.ToString()).ToList());
        }

        private void SetDenovoQualityScores(SampleList<CanvasSegmentsSet> canvasSegmentsSet, SampleList<PedigreeMemberInfo> samplesInfo, List<SampleId> parentIDs, List<SampleId> offspringIDs,
            int canvasSegmentIndex, CopyNumbersLikelihood copyNumberLikelihoods, List<string> names)
        {

            foreach (var probandId in offspringIDs)
            {
                // targeted proband is REF
                if (IsReferenceVariant(canvasSegmentsSet, samplesInfo, probandId, canvasSegmentIndex))
                    continue;
                // common variant
                if (IsCommonCnv(canvasSegmentsSet, samplesInfo, parentIDs, probandId, canvasSegmentIndex))
                    continue;
                // other offsprings are ALT
                if (!offspringIDs.Except(probandId.ToEnumerable()).All(id => IsReferenceVariant(canvasSegmentsSet, samplesInfo, id, canvasSegmentIndex)))
                    continue;
                // not all q-scores are above the threshold
                if (parentIDs.Select(id => !IsPassVariant(canvasSegmentsSet, id, canvasSegmentIndex)).Any() || !IsPassVariant(canvasSegmentsSet, probandId, canvasSegmentIndex))
                    continue;

                double deNovoQualityScore = GetConditionalDeNovoQualityScore(copyNumberLikelihoods, probandId, canvasSegmentsSet, canvasSegmentIndex, names, parentIDs);
                if (Double.IsInfinity(deNovoQualityScore) | deNovoQualityScore > CallerParameters.MaxQscore)
                    deNovoQualityScore = CallerParameters.MaxQscore;
                canvasSegmentsSet[probandId].GetSet()[canvasSegmentIndex].DqScore = deNovoQualityScore;
            }
        }

        private bool IsPassVariant(SampleList<CanvasSegmentsSet> canvasSegmentsSet, SampleId sampleId, int index)
        {
            return canvasSegmentsSet[sampleId].GetSet()[index].QScore > QualityFilterThreshold;
        }

        private bool IsCommonCnv(SampleList<CanvasSegmentsSet> canvasSegmentsSet, SampleList<PedigreeMemberInfo> samplesInfo, List<SampleId> parentIDs, SampleId probandId, int index)
        {
            int parent1CopyNumber = GetCnState(canvasSegmentsSet, parentIDs.First(), index, CallerParameters.MaximumCopyNumber);
            int parent2CopyNumber = GetCnState(canvasSegmentsSet, parentIDs.Last(), index, CallerParameters.MaximumCopyNumber);
            int probandCopyNumber = GetCnState(canvasSegmentsSet, probandId, index, CallerParameters.MaximumCopyNumber);
            var parent1Genotypes = GenerateCnAlleles(parent1CopyNumber);
            var parent2Genotypes = GenerateCnAlleles(parent2CopyNumber);
            var probandGenotypes = GenerateCnAlleles(probandCopyNumber);
            var parent1Segment = canvasSegmentsSet[parentIDs.First()].GetSet()[index];
            var parent2Segment = canvasSegmentsSet[parentIDs.Last()].GetSet()[index];
            var probandSegment = canvasSegmentsSet[probandId].GetSet()[index];
            int parent1Ploidy = samplesInfo[probandId].GetPloidy(parent1Segment);
            int parent2Ploidy = samplesInfo[probandId].GetPloidy(parent2Segment);
            int probandPloidy = samplesInfo[probandId].GetPloidy(probandSegment);
            bool isCommoCnv = parent1Genotypes.Intersect(probandGenotypes).Any() && parent1Ploidy == probandPloidy ||
                              parent2Genotypes.Intersect(probandGenotypes).Any() && parent2Ploidy == probandPloidy;
            return isCommoCnv;
        }

        private bool IsReferenceVariant(SampleList<CanvasSegmentsSet> canvasSegmentsSet, SampleList<PedigreeMemberInfo> samplesInfo, SampleId sampleId, int index)
        {
            var segment = canvasSegmentsSet[sampleId].GetSet()[index];
            return GetCnState(canvasSegmentsSet, sampleId, index, CallerParameters.MaximumCopyNumber) == samplesInfo[sampleId].GetPloidy(segment);
        }


        private void EstimateQScoresNoPedigreeInfo(SampleList<CanvasSegmentsSet> canvasSegmentsSet, int canvasSegmentIndex, double[][] copyNumberLikelihoods)
        {
            var cnStates =
                canvasSegmentsSet.SampleData.Select(
                    x => Math.Min(x.GetSet()[canvasSegmentIndex].CopyNumber, CallerParameters.MaximumCopyNumber - 1)).ToList();
            var counter = 0;
            foreach (var canvasSegmentSet in canvasSegmentsSet.SampleData)
            {
                double normalizationConstant = copyNumberLikelihoods[counter].Sum();
                double qscore = -10.0 * Math.Log10((normalizationConstant - copyNumberLikelihoods[counter][cnStates[counter]]) / normalizationConstant);
                if (Double.IsInfinity(qscore) | qscore > CallerParameters.MaxQscore)
                    qscore = CallerParameters.MaxQscore;
                canvasSegmentSet.GetSet()[canvasSegmentIndex].QScore = qscore;
                if (canvasSegmentSet.GetSet()[canvasSegmentIndex].QScore < QualityFilterThreshold)
                    canvasSegmentSet.GetSet()[canvasSegmentIndex].Filter = $"q{QualityFilterThreshold}";
                counter++;
            }
        }


        private static int GetCnState(SampleList<CanvasSegmentsSet> canvasSegmentsSet, SampleId sampleId, int index, int maximumCopyNumber)
        {
            return Math.Min(canvasSegmentsSet[sampleId].GetSet()[index].CopyNumber, maximumCopyNumber - 1);
        }

        public static int AggregateVariantCoverage(ref List<CanvasSegment> segments)
        {
            var variantCoverage = segments.SelectMany(segment => segment.Balleles.TotalCoverage).ToList();
            return variantCoverage.Any() ? Utilities.Median(variantCoverage) : 0;
        }

        public double[][] GetTransitionMatrix(int numCnStates)
        {
            double[][] transitionMatrix = Utilities.MatrixCreate(numCnStates, numCnStates);
            transitionMatrix[0][0] = 1.0;
            for (int gt = 1; gt < numCnStates; gt++)
                transitionMatrix[0][gt] = 0;

            for (int cn = 1; cn < numCnStates; cn++)
            {
                var gtLikelihood = new Poisson(Math.Max(cn / 2.0, 0.1));
                for (int gt = 0; gt < numCnStates; gt++)
                    transitionMatrix[cn][gt] = gtLikelihood.Probability(gt);
            }
            return transitionMatrix;
        }

        /// <summary>
        /// Identify variant with the highest likelihood at a given setPosition and assign relevant scores
        /// </summary>
        /// <param name="samples"></param>
        /// <param name="setPosition"></param>
        /// <param name="copyNumbers"></param>
        /// <param name="genotypes"></param>
        public void CallVariant(SampleList<CanvasSegmentsSet> canvasSegmentsSet, SampleList<PedigreeMemberInfo> pedigreeMembersInfo, SampleList<CopyNumberModel> model,
             List<List<int>> copyNumbers, Dictionary<int, List<Genotype>> genotypes)
        {
            SegmentsSet segmentSet;

            if (canvasSegmentsSet.SampleData.First().SetA == null)
                segmentSet = SegmentsSet.SetB;
            else if (canvasSegmentsSet.SampleData.First().SetB == null)
                segmentSet = SegmentsSet.SetA;
            else
                segmentSet = GetSegmentSetLikelihoodNoPedigreeInfo(canvasSegmentsSet, pedigreeMembersInfo, model, copyNumbers, SegmentsSet.SetA) >
                              GetSegmentSetLikelihoodNoPedigreeInfo(canvasSegmentsSet, pedigreeMembersInfo, model, copyNumbers, SegmentsSet.SetB) ?
                              SegmentsSet.SetA : SegmentsSet.SetB;

            canvasSegmentsSet.SampleIds.ForEach(id => canvasSegmentsSet[id].SetSet(segmentSet));
            int nSegments = canvasSegmentsSet.SampleData.First().GetSet().Count;

            for (var canvasSegmentIndex = 0; canvasSegmentIndex < nSegments; canvasSegmentIndex++)
            {
                var ll = AssignCopyNumberNoPedigreeInfo(canvasSegmentsSet, pedigreeMembersInfo, model, canvasSegmentIndex, copyNumbers);
                EstimateQScoresNoPedigreeInfo(canvasSegmentsSet, canvasSegmentIndex, ll);
                AssignMccNoPedigreeInfo(canvasSegmentsSet, pedigreeMembersInfo, model, canvasSegmentIndex, genotypes);
            }
        }

        private double GetSegmentSetLikelihoodNoPedigreeInfo(SampleList<CanvasSegmentsSet> canvasSegmentsSet, SampleList<PedigreeMemberInfo> samplesInfo,
            SampleList<CopyNumberModel> copyNumberModel, List<List<int>> copyNumberCombination, SegmentsSet segmentsSet)
        {
            double segmentSetLikelihood = 0;
            foreach (var sampleId in canvasSegmentsSet.SampleIds)
                canvasSegmentsSet[sampleId].SetSet(segmentsSet);
            int nSegments = canvasSegmentsSet.SampleData.First().GetSet().Count;
            for (var canvasSegmentIndex = 0; canvasSegmentIndex < nSegments; canvasSegmentIndex++)
            {
                segmentSetLikelihood += Utilities.MaxValue(AssignCopyNumberNoPedigreeInfo(canvasSegmentsSet, samplesInfo, copyNumberModel,
                    canvasSegmentIndex, copyNumberCombination));
            }
            segmentSetLikelihood /= nSegments;

            return segmentSetLikelihood;
        }

        /// <summary>
        /// Identify variant with the highest likelihood at a given setPosition and assign relevant scores
        /// </summary>
        /// <param name="pedigreeMembers"></param>
        /// <param name="parentIDs"></param>
        /// <param name="children"></param>
        /// <param name="setPosition"></param>
        /// <param name="transitionMatrix"></param>
        /// <param name="offspringsGenotypes"></param>
        /// <param name="genotypes"></param>
        public void CallVariantInPedigree(SampleList<CanvasSegmentsSet> canvasSegmentsSet, SampleList<PedigreeMemberInfo> pedigreeMembersInfo, SampleList<CopyNumberModel> model,
                        List<SampleId> parentIDs, List<SampleId> offspringIDs, double[][] transitionMatrix, List<List<Genotype>> offspringsGenotypes,
                        Dictionary<int, List<Genotype>> genotypes)
        {
            SegmentsSet segmentsSet;

            if (canvasSegmentsSet[parentIDs.First()].SetA == null)
                segmentsSet = SegmentsSet.SetB;
            else if (canvasSegmentsSet[parentIDs.First()].SetB == null)
                segmentsSet = SegmentsSet.SetA;
            else
                segmentsSet = GetSegmentSetLikelihood(canvasSegmentsSet, pedigreeMembersInfo, model, parentIDs,
                                  offspringIDs, SegmentsSet.SetA, transitionMatrix, offspringsGenotypes) >
                              GetSegmentSetLikelihood(canvasSegmentsSet, pedigreeMembersInfo, model, parentIDs,
                                  offspringIDs, SegmentsSet.SetB, transitionMatrix, offspringsGenotypes) ?
                              SegmentsSet.SetA : SegmentsSet.SetB;

            canvasSegmentsSet.SampleIds.ForEach(id => canvasSegmentsSet[id].SetSet(segmentsSet));
            int nSegments = canvasSegmentsSet[parentIDs.First()].GetSet().Count;

            for (var canvasSegmentIndex = 0; canvasSegmentIndex < nSegments; canvasSegmentIndex++)
            {
                var copyNumbersLikelihood = AssignCopyNumberWithPedigreeInfo(canvasSegmentsSet, pedigreeMembersInfo, model, parentIDs,
                    offspringIDs, canvasSegmentIndex, transitionMatrix, offspringsGenotypes);

                EstimateQScoresWithPedigreeInfo(canvasSegmentsSet, pedigreeMembersInfo, parentIDs,
                    offspringIDs, canvasSegmentIndex, copyNumbersLikelihood);

                if (!UseMafInformation(canvasSegmentsSet, canvasSegmentIndex))
                    AssignMccWithPedigreeInfo(canvasSegmentsSet, pedigreeMembersInfo, model, parentIDs,
                        offspringIDs, canvasSegmentIndex, genotypes);
            }
        }

        private double GetSegmentSetLikelihood(SampleList<CanvasSegmentsSet> canvasSegmentsSet,
            SampleList<PedigreeMemberInfo> pedigreeMembersInfo, SampleList<CopyNumberModel> model, List<SampleId> parentIDs,
            List<SampleId> offspringIDs, SegmentsSet segmentsSet, double[][] transitionMatrix,
            List<List<Genotype>> offspringsGenotypes)
        {
            double segmentSetLikelihood = 0;
            foreach (var sampleId in canvasSegmentsSet.SampleIds)
                canvasSegmentsSet[sampleId].SetSet(segmentsSet);

            int nSegments = canvasSegmentsSet.First().Value.GetSet().Count;
            for (var canvasSegmentIndex = 0; canvasSegmentIndex < nSegments; canvasSegmentIndex++)
                segmentSetLikelihood += AssignCopyNumberWithPedigreeInfo(canvasSegmentsSet, pedigreeMembersInfo, model, parentIDs,
                    offspringIDs, canvasSegmentIndex, transitionMatrix, offspringsGenotypes).MaximalLikelihood;

            segmentSetLikelihood /= nSegments;
            return segmentSetLikelihood;
        }


        /// <summary>
        /// Calculates maximal likelihood for copy numbers. Updated CanvasSegment CopyNumber only. 
        /// </summary>
        /// <param name="segmentsSet"></param>
        /// <param name="samplesInfo"></param>
        /// <param name="model"></param>
        /// <param name="parentIDs"></param>
        /// <param name="offspringIDs"></param>
        /// <param name="setPosition"></param>
        /// <param name="segmentPosition"></param>
        /// <param name="segmentsSet"></param>
        /// <param name="transitionMatrix"></param>
        /// <param name="offspringsGenotypes"></param>
        public CopyNumbersLikelihood AssignCopyNumberWithPedigreeInfo(SampleList<CanvasSegmentsSet> segmentsSet,
            SampleList<PedigreeMemberInfo> samplesInfo, SampleList<CopyNumberModel> model, List<SampleId> parentIDs,
            List<SampleId> offspringIDs, int canvasSegmentIndex, double[][] transitionMatrix,
                List<List<Genotype>> offspringsGenotypes)
        {
            int nCopies = CallerParameters.MaximumCopyNumber;
            var names = parentIDs.Union(offspringIDs).Select(x => x.ToString()).ToList();
            var density = new CopyNumbersLikelihood(nCopies, names);
            InitializeCn(segmentsSet, canvasSegmentIndex);
            density.MaximalLikelihood = 0;
            var coverages = parentIDs.Select(id => Math.Min(segmentsSet[id].GetSet()[canvasSegmentIndex].MedianCount,
                samplesInfo[id].MeanCoverage * 3.0)).ToList();
            var parent1Likelihood = model[parentIDs.First()].GetCnLikelihood(coverages.First());
            var parent2Likelihood = model[parentIDs.Last()].GetCnLikelihood(coverages.Last());

            if (parent1Likelihood.Count != parent2Likelihood.Count)
                throw new ArgumentException("Both parentIDs should have the same number of CN states");

            for (int cn1 = 0; cn1 < nCopies; cn1++)
            {
                for (int cn2 = 0; cn2 < nCopies; cn2++)
                {
                    foreach (var offspringGtStates in offspringsGenotypes)
                    {
                        double currentLikelihood = parent1Likelihood[cn1] * parent2Likelihood[cn2];
                        for (var counter = 0; counter < offspringIDs.Count; counter++)
                        {
                            var child = offspringIDs[counter];
                            int modelIndex = Math.Min(offspringGtStates[counter].CountsA + offspringGtStates[counter].CountsB,
                                    CallerParameters.MaximumCopyNumber - 1);
                            double coverage = Math.Min(segmentsSet[child].GetSet()[canvasSegmentIndex].MedianCount, samplesInfo[child].MeanCoverage * 3.0);

                            currentLikelihood *= transitionMatrix[cn1][offspringGtStates[counter].CountsA] *
                                                 transitionMatrix[cn2][offspringGtStates[counter].CountsB] *
                                                 model[child].GetCnLikelihood(coverage)[modelIndex];
                        }
                        int[] copyNumberIndices = { cn1, cn2 };
                        var modelsIndex = copyNumberIndices.Concat(offspringGtStates.Select(x => x.CountsA + x.CountsB)).ToArray();
                        density.SetJointProbability(Math.Max(currentLikelihood, density.GetJointProbability(modelsIndex)), modelsIndex);

                        currentLikelihood = Double.IsNaN(currentLikelihood) || Double.IsInfinity(currentLikelihood)
                            ? 0
                            : currentLikelihood;

                        if (currentLikelihood > density.MaximalLikelihood)
                        {
                            density.MaximalLikelihood = currentLikelihood;
                            segmentsSet[parentIDs.First()].GetSet()[canvasSegmentIndex].CopyNumber = cn1;
                            segmentsSet[parentIDs.Last()].GetSet()[canvasSegmentIndex].CopyNumber = cn2;
                            for (int counter = 0; counter < offspringIDs.Count; counter++)
                            {
                                segmentsSet[offspringIDs[counter]].GetSet()[canvasSegmentIndex].CopyNumber =
                                    offspringGtStates[counter].CountsA + offspringGtStates[counter].CountsB;
                            }
                        }
                    }
                }
            }
            return density;
        }

        /// <summary>
        /// Calculates maximal likelihood for genotypes given a copy number call. Updated MajorChromosomeCount.
        /// </summary>
        /// <param name="parentIDs"></param>
        /// <param name="children"></param>
        /// <param name="segmentPosition"></param>
        /// <param name="genotypes"></param>
        public void AssignMccWithPedigreeInfo(SampleList<CanvasSegmentsSet> canvasSegmentsSet,
            SampleList<PedigreeMemberInfo> samplesInfo, SampleList<CopyNumberModel> model, List<SampleId> parentIDs,
            List<SampleId> offspringIDs, int index, Dictionary<int, List<Genotype>> genotypes)
        {
            double maximalLikelihood = Double.MinValue;
            int parent1CopyNumber = canvasSegmentsSet[parentIDs.First()].GetSet()[index].CopyNumber;
            int parent2CopyNumber = canvasSegmentsSet[parentIDs.Last()].GetSet()[index].CopyNumber;

            foreach (var parent1GtStates in genotypes[parent1CopyNumber])
            {
                foreach (var parent2GtStates in genotypes[parent2CopyNumber])
                {
                    var bestChildGtStates = new List<Genotype>();
                    double currentLikelihood = 1;
                    foreach (SampleId child in offspringIDs)
                    {
                        int childCopyNumber = canvasSegmentsSet[child].GetSet()[index].CopyNumber;
                        bool isInheritedCnv = !canvasSegmentsSet[child].GetSet()[index].DqScore.HasValue;
                        double bestLikelihood = Double.MinValue;
                        Genotype bestGtState = null;
                        bestLikelihood = GetProbandLikelihood(model[child], index, genotypes, childCopyNumber,
                            parent1GtStates, parent2GtStates, isInheritedCnv, canvasSegmentsSet[child], bestLikelihood, ref bestGtState);
                        bestChildGtStates.Add(bestGtState);
                        currentLikelihood *= bestLikelihood;
                    }
                    currentLikelihood *= GetCurrentGtLikelihood(model[parentIDs.First()], canvasSegmentsSet[parentIDs.First()], index, parent1GtStates) *
                                         GetCurrentGtLikelihood(model[parentIDs.Last()], canvasSegmentsSet[parentIDs.Last()], index, parent2GtStates);

                    currentLikelihood = Double.IsNaN(currentLikelihood) || Double.IsInfinity(currentLikelihood)
                        ? 0
                        : currentLikelihood;

                    if (currentLikelihood > maximalLikelihood)
                    {
                        maximalLikelihood = currentLikelihood;
                        AssignMCC(canvasSegmentsSet[parentIDs.First()], model[parentIDs.First()], index, genotypes, parent1GtStates, parent1CopyNumber);
                        AssignMCC(canvasSegmentsSet[parentIDs.Last()], model[parentIDs.Last()], index, genotypes, parent2GtStates, parent2CopyNumber);
                        var counter = 0;
                        foreach (SampleId child in offspringIDs)
                        {
                            if (bestChildGtStates[counter] == null) continue;
                            int childCopyNumber = canvasSegmentsSet[child].GetSet()[index].CopyNumber;
                            AssignMCC(canvasSegmentsSet[child], model[child], index, genotypes, bestChildGtStates[counter], childCopyNumber);
                            counter++;
                        }
                    }
                }
            }
        }

        private double GetProbandLikelihood(CopyNumberModel copyNumberModel, int index, Dictionary<int, List<Genotype>> genotypes,
            int childCopyNumber, Genotype parent1GtStates, Genotype parent2GtStates, bool isInheritedCnv, CanvasSegmentsSet canvasSegmentsSet,
            double bestLikelihood, ref Genotype bestGtState)
        {
            foreach (var childGtState in genotypes[childCopyNumber])
            {
                double currentChildLikelihood;
                if (IsGtPedigreeConsistent(parent1GtStates, childGtState) &&
                    IsGtPedigreeConsistent(parent2GtStates, childGtState)
                    && isInheritedCnv)
                    currentChildLikelihood = copyNumberModel.GetCurrentGtLikelihood(canvasSegmentsSet.GetSet()[index].Balleles.GetAlleleCounts(), childGtState);
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

        private static void AssignMCC(CanvasSegmentsSet canvasSegmentsSet, CopyNumberModel copyNumberModel, int index,
            Dictionary<int, List<Genotype>> genotypes, Genotype gtStates, int copyNumber)
        {
            const int diploidCopyNumber = 2;
            const int haploidCopyNumber = 1;
            if (copyNumber > diploidCopyNumber)
            {

                canvasSegmentsSet.GetSet()[index].MajorChromosomeCount =
                    Math.Max(gtStates.CountsA, gtStates.CountsB);
                int? selectedGtState = genotypes[copyNumber].IndexOf(gtStates);
                canvasSegmentsSet.GetSet()[index].MajorChromosomeCountScore =
                    copyNumberModel.GetGtLikelihoodScore(canvasSegmentsSet.GetSet()[index].Balleles.GetAlleleCounts(), genotypes[copyNumber], ref selectedGtState);
            }
            else
            {
                canvasSegmentsSet.GetSet()[index].MajorChromosomeCount = copyNumber == diploidCopyNumber
                    ? haploidCopyNumber : copyNumber;
                canvasSegmentsSet.GetSet()[index].MajorChromosomeCountScore = null;
            }
        }

        private static double GetCurrentGtLikelihood(CopyNumberModel copyNumberModel, CanvasSegmentsSet canvasSegmentsSet, int index, Genotype gtStates)
        {
            return copyNumberModel.GetCurrentGtLikelihood(canvasSegmentsSet.GetSet()[index].Balleles.GetAlleleCounts(), gtStates);
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
        /// <param name="canvasSegmentsSet"></param>
        /// <param name="segmentPosition"></param>
        /// <param name="copyNumberCombinations"></param>
        public double[][] AssignCopyNumberNoPedigreeInfo(SampleList<CanvasSegmentsSet> canvasSegmentsSet, SampleList<PedigreeMemberInfo> samplesInfo,
            SampleList<CopyNumberModel> copyNumberModel, int canvasSegmentIndex, List<List<int>> copyNumberCombinations)
        {
            const int defaultCn = 2;
            const double maxCoverageMultiplier = 3.0;

            double maximalLikelihood;
            foreach (var canvasSegmentSet in canvasSegmentsSet.SampleData)
                canvasSegmentSet.GetSet()[canvasSegmentIndex].CopyNumber = defaultCn;
            int nCopies = CallerParameters.MaximumCopyNumber;
            var names = canvasSegmentsSet.SampleIds.Select(x => x.ToString()).ToList();
            var totalLikelihoods = new List<double>();
            foreach (var copyNumberCombination in copyNumberCombinations)
            {
                double totalLikelihood = 0;
                foreach (var sampleId in canvasSegmentsSet.SampleIds)
                {
                    maximalLikelihood = 0;
                    foreach (var copyNumber in copyNumberCombination)
                    {
                        double currentLikelihood =
                            copyNumberModel[sampleId].GetCnLikelihood(
                                Math.Min(canvasSegmentsSet[sampleId].GetSet()[canvasSegmentIndex].MedianCount,
                                    samplesInfo[sampleId].MeanCoverage * maxCoverageMultiplier))[copyNumber];
                        currentLikelihood = Double.IsNaN(currentLikelihood) || Double.IsInfinity(currentLikelihood)
                            ? 0
                            : currentLikelihood;
                        if (currentLikelihood > maximalLikelihood)
                            maximalLikelihood = currentLikelihood;
                    }
                    totalLikelihood += maximalLikelihood;
                }
                totalLikelihoods.Add(totalLikelihood);
            }

            var density = new double[canvasSegmentsSet.Count()][];
            // no need to iterate over multiple genotypes for n (segmentsSet.Count) = 1
            if (canvasSegmentsSet.Count() == 1)
            {
                canvasSegmentsSet.Single().Value.GetSet()[canvasSegmentIndex].CopyNumber =
                    copyNumberCombinations[totalLikelihoods.IndexOf(totalLikelihoods.Max())].First();
                density[0] = totalLikelihoods.ToArray();
                return density;
            }

            var bestCopyNumberCombination = copyNumberCombinations[totalLikelihoods.IndexOf(totalLikelihoods.Max())];
            int counter = 0;
            foreach (var sampleId in canvasSegmentsSet.SampleIds)
            {
                maximalLikelihood = 0;
                density[counter] = new double[nCopies];
                counter++;
                foreach (var copyNumber in bestCopyNumberCombination)
                {
                    double currentLikelihood =
                        copyNumberModel[sampleId].GetCnLikelihood(
                            Math.Min(canvasSegmentsSet[sampleId].GetSet()[canvasSegmentIndex].MedianCount,
                                samplesInfo[sampleId].MeanCoverage * maxCoverageMultiplier))[copyNumber];
                    currentLikelihood = Double.IsNaN(currentLikelihood) || Double.IsInfinity(currentLikelihood)
                        ? 0
                        : currentLikelihood;
                    if (currentLikelihood > maximalLikelihood)
                    {
                        maximalLikelihood = currentLikelihood;
                        canvasSegmentsSet[sampleId].GetSet()[canvasSegmentIndex].CopyNumber = copyNumber;
                        density[names.FindIndex(name => name == sampleId.ToString())][copyNumber] = maximalLikelihood;
                    }
                }
            }
            return density;
        }

        /// <summary>
        /// Calculates maximal likelihood for segments with SNV allele counts given CopyNumber. Updated MajorChromosomeCount.
        /// </summary>
        /// <param name="samples"></param>
        /// <param name="segmentPosition"></param>
        /// <param name="genotypes"></param>       
        public void AssignMccNoPedigreeInfo(SampleList<CanvasSegmentsSet> canvasSegmentsSet, SampleList<PedigreeMemberInfo> pedigreeMembersInfo,
            SampleList<CopyNumberModel> model, int canvasSegmentIndex, Dictionary<int, List<Genotype>> genotypes)
        {
            foreach (var sampleId in canvasSegmentsSet.SampleIds)
            {
                int copyNumber = canvasSegmentsSet[sampleId].GetSet()[canvasSegmentIndex].CopyNumber;
                if (copyNumber > 2)
                {
                    canvasSegmentsSet[sampleId].GetSet()[canvasSegmentIndex].MajorChromosomeCount = copyNumber == 2 ? 1 : copyNumber;
                    return;
                }
                var genotypeset = genotypes[copyNumber];
                int? selectedGtState = null;
                double gqscore = model[sampleId].GetGtLikelihoodScore(canvasSegmentsSet[sampleId].GetSet()[canvasSegmentIndex].Balleles.GetAlleleCounts(),
                    genotypeset, ref selectedGtState);
                canvasSegmentsSet[sampleId].GetSet()[canvasSegmentIndex].MajorChromosomeCountScore = gqscore;
                if (selectedGtState.HasValue)
                    canvasSegmentsSet[sampleId].GetSet()[canvasSegmentIndex].MajorChromosomeCount =
                        Math.Max(genotypeset[selectedGtState.Value].CountsA,
                            genotypeset[selectedGtState.Value].CountsB);
            }
        }


        private static void InitializeCn(SampleList<CanvasSegmentsSet> samples, int canvasSegmentIndex)
        {
            const int defaultCn = 2;
            samples.ForEach(sample => sample.Value.GetSet()[canvasSegmentIndex].CopyNumber = defaultCn);
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
                var permutations = new Combinations<int>(cnStates, currentAlleleNumber);
                var list = permutations.Select(x => x.ToList()).ToList();
                allCombinations.AddRange(list);
            }
            return allCombinations;
        }

        public List<Genotype> GenerateParentalGenotypes(int numCnStates)
        {
            var genotypes = new List<Genotype>();
            for (int cn = 0; cn < numCnStates; cn++)
            {
                for (int gt = 0; gt <= cn; gt++)
                {
                    genotypes.Add(new Genotype(gt, cn - gt));
                }
            }
            return genotypes;
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

        public void GenerateOffspringGenotypes(List<List<Genotype>> offspringGenotypes, List<Genotype> genotypeSet,
            int nOffsprings, List<Genotype> partialGenotypes)
        {
            if (nOffsprings > 0)
            {
                foreach (var genotype in genotypeSet)
                {
                    GenerateOffspringGenotypes(offspringGenotypes, genotypeSet, nOffsprings - 1,
                        partialGenotypes.Concat(new List<Genotype> { genotype }).ToList());
                }
            }
            if (nOffsprings == 0)
            {
                offspringGenotypes.Add(partialGenotypes);
            }
        }

        private static List<SegmentIndexRange> GetParallelIntervals(int nSegments, int nCores)
        {

            var segmentIndexRanges = new List<SegmentIndexRange>();

            int step = nSegments / nCores;
            segmentIndexRanges.Add(new SegmentIndexRange(0, step));
            int cumSum = step + 1;
            while (cumSum + step + 1 < nSegments - 1)
            {
                segmentIndexRanges.Add(new SegmentIndexRange(cumSum, cumSum + step));
                cumSum += step + 1;
            }
            segmentIndexRanges.Add(new SegmentIndexRange(cumSum, nSegments - 1));
            return segmentIndexRanges;
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
            Other, Parent, Proband
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

        public double GetSingleSampleQualityScore(CopyNumbersLikelihood density, int cnState, string sampleName)
        {
            var cnMarginalProbabilities = density.GetMarginalProbability(CallerParameters.MaximumCopyNumber, sampleName);
            double normalizationConstant = cnMarginalProbabilities.Sum();
            double qscore = -10.0 * Math.Log10((normalizationConstant - cnMarginalProbabilities[cnState]) / normalizationConstant);
            if (Double.IsInfinity(qscore) | qscore > CallerParameters.MaxQscore)
                qscore = CallerParameters.MaxQscore;

            return qscore;
        }

        public double GetConditionalDeNovoQualityScore(CopyNumbersLikelihood density, SampleId probandId, SampleList<CanvasSegmentsSet> canvasSegmentsSet,
            int index, List<string> names, List<SampleId> parentIDs)
        {

            var numerator = 0.0;
            var denominator = 0.0;
            const int diploidState = 2;
            int parent1Index = names.IndexOf(parentIDs.First().ToString());
            int parent2Index = names.IndexOf(parentIDs.Last().ToString());
            int probandIndex = names.IndexOf(probandId.ToString());

            var probandMarginalProbabilities = density.GetMarginalProbability(CallerParameters.MaximumCopyNumber, probandId.ToString());
            int probandCopyNumber = GetCnState(canvasSegmentsSet, probandId, index, CallerParameters.MaximumCopyNumber);
            double normalization = probandMarginalProbabilities[probandCopyNumber] + probandMarginalProbabilities[diploidState];
            double probandMarginalAlt = probandMarginalProbabilities[probandCopyNumber] / normalization;

            foreach (var copyNumberIndex in density.Indices.Where(x => x[probandIndex] == probandCopyNumber))
            {
                if (!(density.GetJointProbability(copyNumberIndex.ToArray()) > 0.0))
                    continue;
                double holder = density.GetJointProbability(copyNumberIndex);
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



