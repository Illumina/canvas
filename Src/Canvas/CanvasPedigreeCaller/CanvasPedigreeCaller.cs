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
            string outVcfFile, string ploidyBedPath,
            string referenceFolder, List<string> sampleNames, string commonCNVsbedPath, string pedigreeFile)
        {
            // load files
            // initialize data structures and classes
            int fileCounter = 0;
            var kinships = ReadPedigreeFile(pedigreeFile);
            var pedigreeMembers = new SampleList<PedigreeMember>();
            var pedigreeMembersInfo = new SampleList<PedigreeMemberInfo>();
            var sampleSegments = new SampleList<Segments>();
            var copyNumberModels = new SampleList<CopyNumberModel>();


            foreach (string sampleName in sampleNames)
            {
                var sampleId = new SampleId(sampleName);
                var segment = Segments.ReadSegments(_logger, new FileLocation(segmentFiles[fileCounter]));
                segment.AddAlleles(CanvasIO.ReadFrequenciesWrapper(_logger, new FileLocation(variantFrequencyFiles[fileCounter]), segment.IntervalsByChromosome));
                sampleSegments.Add(sampleId, segment);

                var pedigreeMember = GetPedigreeMember(sampleSegments[sampleId], variantFrequencyFiles[fileCounter], segmentFiles[fileCounter], sampleId,
                    CallerParameters.DefaultReadCountsThreshold, kinships[sampleId], commonCNVsbedPath);
                var pedigreeMemberInfo = PedigreeMemberInfo.GetPedigreeMemberInfo(segment, ploidyBedPath, CallerParameters.NumberOfTrimmedBins, sampleId);
                var copyNumberModel = new CopyNumberModel(CallerParameters.MaximumCopyNumber, pedigreeMemberInfo);
                pedigreeMembers.Add(sampleId, pedigreeMember);
                pedigreeMembersInfo.Add(sampleId, pedigreeMemberInfo);
                copyNumberModels.Add(sampleId, copyNumberModel);
                fileCounter++;
            }

            int numberOfSegments = sampleSegments.SampleData.First().AllSegments.Count;
            var segmentIntervals = GetParallelIntervals(numberOfSegments, Math.Min(Environment.ProcessorCount, CallerParameters.MaxCoreNumber));
            var parentsIds = kinships.Where(kin => kin.Value.Equals(PedigreeMember.Kinship.Parent)).Select(kin => kin.Key).ToList();
            var offspringsIds = kinships.Where(kin => kin.Value.Equals(PedigreeMember.Kinship.Proband)).Select(kin => kin.Key).ToList();
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
                segmentIntervals,
                interval =>
                {
                    Console.WriteLine($"{DateTime.Now} Launching SPW task for segment {interval.Start} - {interval.End}");
                    for (int segmentIndex = interval.Start; segmentIndex <= interval.End; segmentIndex++)
                    {
                        CallVariantInPedigree(pedigreeMembers, pedigreeMembersInfo, copyNumberModels, parentsIds,
                            offspringsIds, segmentIndex, transitionMatrix, offspringsGenotypes, genotypes);
                    }
                    Console.WriteLine($"{DateTime.Now} Finished SPW task for segment {interval.Start} - {interval.End}");
                });

            var variantCalledSegments = new SampleList<List<CanvasSegment>>();
            foreach (var key in pedigreeMembers.SampleIds)
                variantCalledSegments.Add(key, pedigreeMembers[key].GetCanvasSegments());

            var mergedVariantCalledSegments = MergeSegments(variantCalledSegments, CallerParameters.MinimumCallSize);
            var outputFolder = new FileLocation(outVcfFile).Directory;
            foreach (var sampleId in offspringsIds.Union(parentsIds))
            {
                var coverageOutputPath = SingleSampleCallset.GetCoverageAndVariantFrequencyOutput(outputFolder,
                    sampleId.ToString());
                CanvasSegment.WriteCoveragePlotData(mergedVariantCalledSegments[sampleId], pedigreeMembersInfo[sampleId].MeanCoverage,
                    pedigreeMembersInfo[sampleId].Ploidy, coverageOutputPath.FullName, referenceFolder);
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
        /// <param name="pedigreeMembers"></param>
        /// <param name="segmentIndex"></param>
        /// <returns></returns>
        private bool UseMafInformation(SampleList<PedigreeMember> pedigreeMembers, int segmentSetIndex, int segmentIndex, SegmentsSet segmentsSet)
        {
            var alleles = pedigreeMembers.SampleData.Select(x => x.SegmentSets[segmentSetIndex].GetSet(segmentsSet)[segmentIndex].Balleles?.TotalCoverage);
            var alleleCounts = alleles.Select(allele => allele?.Count ?? 0).ToList();
            bool lowAlleleCounts = alleleCounts.Select(x => x < CallerParameters.DefaultReadCountsThreshold).Any(c => c == true);
            var coverageCounts = pedigreeMembers.SampleData.Select(x => x.SegmentSets[segmentSetIndex].GetSet(segmentsSet)[segmentIndex].MedianCount).ToList();
            var isSkewedHetHomRatio = false;
            double alleleDensity = pedigreeMembers.SampleData.First().SegmentSets[segmentSetIndex].GetSet(segmentsSet)[segmentIndex].Length /
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
            const PedigreeMember.Kinship kinship = PedigreeMember.Kinship.Other;
            var pedigreeMembers = new SampleList<PedigreeMember>();
            var pedigreeMembersInfo = new SampleList<PedigreeMemberInfo>();
            var sampleSegments = new SampleList<Segments>();
            var copyNumberModels = new SampleList<CopyNumberModel>();

            foreach (string sampleName in sampleNames)
            {
                var sampleId = new SampleId(sampleName);
                var segment = Segments.ReadSegments(_logger, new FileLocation(segmentFiles[fileCounter]));
                segment.AddAlleles(CanvasIO.ReadFrequenciesWrapper(_logger, new FileLocation(variantFrequencyFiles[fileCounter]), segment.IntervalsByChromosome));
                sampleSegments.Add(sampleId, segment);
                var pedigreeMember = GetPedigreeMember(sampleSegments[sampleId], variantFrequencyFiles[fileCounter], segmentFiles[fileCounter], sampleId,
                    CallerParameters.DefaultReadCountsThreshold, kinship, commonCNVsbedPath);
                var pedigreeMemberInfo = PedigreeMemberInfo.GetPedigreeMemberInfo(segment, ploidyBedPath, CallerParameters.NumberOfTrimmedBins, sampleId);
                var copyNumberModel = new CopyNumberModel(CallerParameters.MaximumCopyNumber, pedigreeMemberInfo);
                pedigreeMembers.Add(sampleId, pedigreeMember);
                pedigreeMembersInfo.Add(sampleId, pedigreeMemberInfo);
                copyNumberModels.Add(sampleId, copyNumberModel);
                fileCounter++;
            }

            int numberOfSegments = sampleSegments.SampleData.First().AllSegments.Count;
            var segmentIntervals = GetParallelIntervals(numberOfSegments, Math.Min(Environment.ProcessorCount, CallerParameters.MaxCoreNumber));
            var genotypes = GenerateGenotypeCombinations(CallerParameters.MaximumCopyNumber);
            int maxAlleleNumber = Math.Min(CallerParameters.MaxAlleleNumber, pedigreeMembers.Count());
            var copyNumberCombinations = GenerateCopyNumberCombinations(CallerParameters.MaximumCopyNumber,
                maxAlleleNumber);

            Parallel.ForEach(
                segmentIntervals,
                interval =>
                {
                    Console.WriteLine($"{DateTime.Now} Launching SPW task for segment {interval.Start} - {interval.End}");
                    for (int segmentIndex = interval.Start; segmentIndex <= interval.End; segmentIndex++)
                    {
                        CallVariant(pedigreeMembers, pedigreeMembersInfo, copyNumberModels, segmentIndex, copyNumberCombinations, genotypes);
                    }
                    Console.WriteLine($"{DateTime.Now} Finished SPW task for segment {interval.Start} - {interval.End}");
                });


            var variantCalledSegments = new SampleList<List<CanvasSegment>>();
            foreach (var key in pedigreeMembers.SampleIds)
                variantCalledSegments.Add(key, pedigreeMembers[key].GetCanvasSegments());

            var mergedVariantCalledSegments = MergeSegments(variantCalledSegments, CallerParameters.MinimumCallSize);
            var outputFolder = new FileLocation(outVcfFile).Directory;
            foreach (string sampleName in sampleNames)
            {
                var sampleId = new SampleId(sampleName);
                var coverageOutputPath = SingleSampleCallset.GetCoverageAndVariantFrequencyOutput(outputFolder, sampleName);
                CanvasSegment.WriteCoveragePlotData(mergedVariantCalledSegments[sampleId], pedigreeMembersInfo[sampleId].MeanCoverage,
                    pedigreeMembersInfo[sampleId].Ploidy, coverageOutputPath.FullName, referenceFolder);
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

        private PedigreeMember GetPedigreeMember(Segments segments, string variantFrequencyFile, string segmentFile, SampleId sampleId, int defaultAlleleCountThreshold, PedigreeMember.Kinship kinship, string commonCNVsbedPath)
        {
            var pedigreeMember = new PedigreeMember(sampleId, kinship);

            if (commonCNVsbedPath != null)
            {
                var segmentsByChromosome = CanvasSegment.GetSegmentsByChromosome(segments.AllSegments);
                var coverage = CanvasSegment.ReadBEDInput(segmentFile);
                var commonRegions = Utilities.LoadBedFile(commonCNVsbedPath);
                Utilities.SortAndOverlapCheck(commonRegions, commonCNVsbedPath);
                if (IdenticalChromosomeNames(commonRegions, coverage) == 0)
                    throw new ArgumentException(
                        $"Chromosome names in a common CNVs bed file {commonCNVsbedPath} does not match " +
                        $"chromosomes in {segmentFile}");

                var segmentIntervalsByChromosome = new Dictionary<string, List<BedInterval>>();
                Parallel.ForEach(commonRegions.Keys, chr => segmentIntervalsByChromosome[chr] =
                CanvasSegment.RemapCommonRegions(commonRegions[chr], coverage.StartByChr[chr], coverage.EndByChr[chr]));
                var allelesByChromosomeCommonSegs = CanvasIO.ReadFrequenciesWrapper(_logger, new FileLocation(variantFrequencyFile), segmentIntervalsByChromosome);
                var segmentsSetByChromosome = new ConcurrentDictionary<string, List<CanvasSegmentsSet>>();
                Parallel.ForEach(
                    segmentsByChromosome.Keys,
                    chr =>
                    {
                        Console.WriteLine($"SegmentsFromCommonCnvs for {chr} started");

                        if (commonRegions.Keys.Any(chromosome => chromosome == chr))
                        {
                            Console.WriteLine($"CreateSegmentsFromCommonCnvs for {chr} ");
                            var commonCnvCanvasSegments = CanvasSegment.CreateSegmentsFromCommonCnvs(coverage, chr,
                                segmentIntervalsByChromosome[chr]);
                            for (var index = 0; index < commonCnvCanvasSegments.Count; index++)
                            {
                                commonCnvCanvasSegments[index].Balleles.Add(allelesByChromosomeCommonSegs[chr][index]);
                            }
                            segmentsSetByChromosome[chr] = CanvasSegment.MergeCommonCnvSegments(segmentsByChromosome[chr], commonCnvCanvasSegments, chr, defaultAlleleCountThreshold) ??
                            segmentsByChromosome[chr].Select(segment => new CanvasSegmentsSet(setA: new List<CanvasSegment> { segment }, setB: null)).ToList();

                            Console.WriteLine($"SegmentsFromCommonCnvs for {chr} count {segmentsByChromosome[chr].Count}");

                            Console.WriteLine($"SegmentsFromCommonCnvs for {chr} returned");
                        }
                        else
                        {
                            segmentsSetByChromosome[chr] = segmentsByChromosome[chr].Select(segment =>
                            new CanvasSegmentsSet(setA: new List<CanvasSegment> { segment }, setB: null)).ToList();
                            Console.WriteLine($"SegmentsFromCommonCnvs for {chr} count {segmentsByChromosome[chr].Count}");
                            Console.WriteLine($"SegmentsFromCommonCnvs for {chr} returned");
                        }
                        Console.WriteLine($"SegmentsFromCommonCnvs for {chr} finished");

                    });
                Console.WriteLine($"SegmentSets for {sampleId} ");
                pedigreeMember.SegmentSets.AddRange(segmentsSetByChromosome.OrderBy(i => i.Key).Select(x => x.Value).SelectMany(x => x).ToList());
            }
            else
            {
                pedigreeMember.SegmentSets =
                    segments.AllSegments.Select(
                            segment =>
                                new CanvasSegmentsSet(setA: new List<CanvasSegment> { segment }, setB: null))
                        .ToList();
            }

            return pedigreeMember;
        }

        private static int IdenticalChromosomeNames(Dictionary<string, List<SampleGenomicBin>> commonRegions,
            CoverageInfo coverage)
        {
            var chromsomes = new HashSet<string>(coverage.CoverageByChr.Keys);
            return commonRegions.Keys.Count(chromosome => chromsomes.Contains(chromosome));
        }

        private void EstimateQScoresWithPedigreeInfo(SampleList<PedigreeMember> pedigreeMembers, SampleList<PedigreeMemberInfo> pedigreeMembersInfo, SampleList<CopyNumberModel> model,
            List<SampleId> parentIDs, List<SampleId> offspringIDs, CanvasSegmentIndex index, CopyNumbersLikelihood copyNumberLikelihoods)
        {
            var cnStates = GetCnStates(pedigreeMembers, parentIDs, offspringIDs, index, CallerParameters.MaximumCopyNumber);
            var names = offspringIDs.Union(parentIDs).Select(x => x.ToString()).ToList();

            var singleSampleQualityScores = GetSingleSampleQualityScores(copyNumberLikelihoods, cnStates, names);
            var counter = 0;
            foreach (var sample in pedigreeMembers)
            {
                sample.Value.GetCanvasSegment(index).QScore = singleSampleQualityScores[counter];
                if (sample.Value.GetCanvasSegment(index).QScore < QualityFilterThreshold)
                    sample.Value.GetCanvasSegment(index).Filter = $"q{QualityFilterThreshold}";
                counter++;
            }
            SetDenovoQualityScores(pedigreeMembers, pedigreeMembersInfo, parentIDs, offspringIDs, index,
                copyNumberLikelihoods, names, cnStates, singleSampleQualityScores);
        }

        private void SetDenovoQualityScores(SampleList<PedigreeMember> samples, SampleList<PedigreeMemberInfo> samplesInfo, List<SampleId> parentIDs, List<SampleId> offspringIDs,
            CanvasSegmentIndex index, CopyNumbersLikelihood copyNumberLikelihoods, List<string> names, List<int> cnStates,
            List<double> singleSampleQualityScores)
        {
            int parent1Index = names.IndexOf(parentIDs.First().ToString());
            int parent2Index = names.IndexOf(parentIDs.Last().ToString());

            foreach (var probandId in offspringIDs)
            {
                int probandIndex = names.IndexOf(probandId.ToString());
                var remainingProbandIndex = offspringIDs.Except(probandId.ToEnumerable()).Select(x => names.IndexOf(x.ToString())).ToList();
                var probandSegment = samples[probandId].GetCanvasSegment(index);
                var probandPloidy = samplesInfo[probandId].GetPloidy(probandSegment);
                if (cnStates[probandIndex] != probandPloidy &&
                    // targeted proband is ALT
                    (ParentsRefCheck(samples, samplesInfo, parentIDs, index, cnStates, parent1Index, parent2Index) ||
                     // either parent are REF or 
                     IsNotCommonCnv(samples, samplesInfo, parentIDs, cnStates, parent1Index, parent2Index, probandIndex, index, probandPloidy)) &&
                    // or a common variant 
                    remainingProbandIndex.Zip(offspringIDs.Except(probandId.ToEnumerable())).All(i =>
                            cnStates[i.Item1] == samplesInfo[i.Item2].GetPloidy(samples[i.Item2].GetCanvasSegment(index)) ||
                            IsNotCommonCnv(samples, samplesInfo, parentIDs, cnStates, parent1Index, parent2Index, probandIndex, index, probandPloidy)) &&
                    // and other probands are REF or common variant 
                    singleSampleQualityScores[probandIndex] > QualityFilterThreshold &&
                    singleSampleQualityScores[parent1Index] > QualityFilterThreshold &&
                    singleSampleQualityScores[parent2Index] > QualityFilterThreshold)
                // and all q-scores are above the threshold
                {
                    double deNovoQualityScore = GetConditionalDeNovoQualityScore(copyNumberLikelihoods, probandIndex,
                        cnStates[probandIndex], names[probandIndex], parent1Index, parent2Index,
                        remainingProbandIndex.ToList());
                    if (Double.IsInfinity(deNovoQualityScore) | deNovoQualityScore > CallerParameters.MaxQscore)
                        deNovoQualityScore = CallerParameters.MaxQscore;
                    samples[probandId].GetCanvasSegment(index).DqScore =
                        deNovoQualityScore;
                }
            }
        }

        private static bool IsNotCommonCnv(SampleList<PedigreeMember> samples, SampleList<PedigreeMemberInfo> samplesInfo, List<SampleId> parentIDs, List<int> cnStates,
            int parent1Index, int parent2Index, int probandIndex, CanvasSegmentIndex index, int probandPloidy)
        {
            var parent1Genotypes = GenerateCnAlleles(cnStates[parent1Index]);
            var parent2Genotypes = GenerateCnAlleles(cnStates[parent2Index]);
            var probandGenotypes = GenerateCnAlleles(cnStates[probandIndex]);
            var segment1 = samples[parentIDs.First()].GetCanvasSegment(index);
            var segment2 = samples[parentIDs.Last()].GetCanvasSegment(index);
            bool isCommoCnv = (parent1Genotypes.Intersect(probandGenotypes).Any() &&
                               samplesInfo[parentIDs.First()].GetPloidy(segment1) == probandPloidy) ||
                              (parent2Genotypes.Intersect(probandGenotypes).Any() &&
                               samplesInfo[parentIDs.Last()].GetPloidy(segment2) == probandPloidy);
            return !isCommoCnv;
        }

        private static bool ParentsRefCheck(SampleList<PedigreeMember> samples, SampleList<PedigreeMemberInfo> samplesInfo, List<SampleId> parentIDs, CanvasSegmentIndex index,
            List<int> cnStates, int parent1Index, int parent2Index)
        {
            var segment1 = samples[parentIDs.First()].GetCanvasSegment(index);
            var segment2 = samples[parentIDs.Last()].GetCanvasSegment(index);
            return cnStates[parent1Index] == samplesInfo[parentIDs.First()].GetPloidy(segment1) &&
                   cnStates[parent2Index] == samplesInfo[parentIDs.Last()].GetPloidy(segment2);
        }

        private void EstimateQScoresNoPedigreeInfo(SampleList<PedigreeMember> pedigreeMembers, CanvasSegmentIndex index, double[][] copyNumberLikelihoods)
        {
            var cnStates =
                pedigreeMembers.SampleData.Select(
                    x => Math.Min(x.GetCanvasSegment(index).CopyNumber, CallerParameters.MaximumCopyNumber - 1)).ToList();
            int counter = 0;
            foreach (PedigreeMember sample in pedigreeMembers.SampleData)
            {
                double normalizationConstant = copyNumberLikelihoods[counter].Sum();
                double qscore = -10.0 * Math.Log10((normalizationConstant - copyNumberLikelihoods[counter][cnStates[counter]]) / normalizationConstant);
                if (Double.IsInfinity(qscore) | qscore > CallerParameters.MaxQscore)
                    qscore = CallerParameters.MaxQscore;
                sample.GetCanvasSegment(index).QScore = qscore;
                if (sample.GetCanvasSegment(index).QScore < QualityFilterThreshold)
                    sample.GetCanvasSegment(index).Filter = $"q{QualityFilterThreshold}";
                counter++;
            }
        }


        private static List<int> GetCnStates(SampleList<PedigreeMember> samples, List<SampleId> parentIDs, List<SampleId> offspringIDs,
            CanvasSegmentIndex index, int maximumCopyNumber)
        {
            return parentIDs.Select(id => Math.Min(samples[id].GetCanvasSegment(index).CopyNumber, maximumCopyNumber - 1))
                    .Concat(offspringIDs.Select(id => Math.Min(samples[id].GetCanvasSegment(index).CopyNumber, maximumCopyNumber - 1))).ToList();
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
        public void CallVariant(SampleList<PedigreeMember> pedigreeMembers, SampleList<PedigreeMemberInfo> pedigreeMembersInfo, SampleList<CopyNumberModel> model,
            int setPosition, List<List<int>> copyNumbers, Dictionary<int, List<Genotype>> genotypes)
        {
            SegmentsSet segmentsSet;

            if (pedigreeMembers.SampleData.First().SegmentSets[setPosition].SetA == null)
                segmentsSet = SegmentsSet.SetB;
            else if (pedigreeMembers.SampleData.First().SegmentSets[setPosition].SetB == null)
                segmentsSet = SegmentsSet.SetA;
            else
                segmentsSet = GetSegmentSetLikelihoodNoPedigreeInfo(pedigreeMembers, pedigreeMembersInfo, model, setPosition, copyNumbers, SegmentsSet.SetA) >
                              GetSegmentSetLikelihoodNoPedigreeInfo(pedigreeMembers, pedigreeMembersInfo, model, setPosition, copyNumbers, SegmentsSet.SetB) ?
                              SegmentsSet.SetA : SegmentsSet.SetB;

            pedigreeMembers.SampleData.ForEach(sample => sample.SegmentSets[setPosition].SelectedSet = segmentsSet);

            for (var segmentPosition = 0;
                segmentPosition <
                pedigreeMembers.SampleData.First().SegmentSets[setPosition].GetSet(segmentsSet).Count;
                segmentPosition++)
            {
                var canvasSegmentIndex = new CanvasSegmentIndex(setPosition, segmentsSet, segmentPosition);

                var ll = AssignCopyNumberNoPedigreeInfo(pedigreeMembers, pedigreeMembersInfo, model, canvasSegmentIndex, copyNumbers);
                EstimateQScoresNoPedigreeInfo(pedigreeMembers, canvasSegmentIndex, ll);
                AssignMccNoPedigreeInfo(pedigreeMembers, pedigreeMembersInfo, model, canvasSegmentIndex, genotypes);
            }
        }

        private double GetSegmentSetLikelihoodNoPedigreeInfo(SampleList<PedigreeMember> samples, SampleList<PedigreeMemberInfo> samplesInfo,
            SampleList<CopyNumberModel> copyNumberModel, int setPosition,
            List<List<int>> copyNumberCombination, SegmentsSet segmentsSet)
        {
            double segmentSetLikelihood = 0;
            int nSegments = samples.SampleData.First().SegmentSets[setPosition].GetSet(segmentsSet).Count;
            for (var segmentPosition = 0; segmentPosition < nSegments; segmentPosition++)
            {
                var canvasSegmentIndex = new CanvasSegmentIndex(setPosition, segmentsSet, segmentPosition);
                segmentSetLikelihood += Utilities.MaxValue(
                    AssignCopyNumberNoPedigreeInfo(samples, samplesInfo, copyNumberModel, canvasSegmentIndex, copyNumberCombination));
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
        public void CallVariantInPedigree(SampleList<PedigreeMember> pedigreeMembers, SampleList<PedigreeMemberInfo> pedigreeMembersInfo, SampleList<CopyNumberModel> model,
                        List<SampleId> parentIDs, List<SampleId> offspringIDs, int setPosition, double[][] transitionMatrix, List<List<Genotype>> offspringsGenotypes,
                        Dictionary<int, List<Genotype>> genotypes)
        {
            SegmentsSet segmentsSet;

            if (pedigreeMembers[parentIDs.First()].SegmentSets[setPosition].SetA == null)
                segmentsSet = SegmentsSet.SetB;
            else if (pedigreeMembers[parentIDs.First()].SegmentSets[setPosition].SetB == null)
                segmentsSet = SegmentsSet.SetA;
            else
                segmentsSet = GetSegmentSetLikelihood(pedigreeMembers, pedigreeMembersInfo, model, parentIDs,
                                  offspringIDs, setPosition, SegmentsSet.SetA, transitionMatrix, offspringsGenotypes) >
                              GetSegmentSetLikelihood(pedigreeMembers, pedigreeMembersInfo, model, parentIDs,
                                  offspringIDs, setPosition, SegmentsSet.SetB, transitionMatrix, offspringsGenotypes) ?
                              SegmentsSet.SetA : SegmentsSet.SetB;

            parentIDs.ForEach(sampleId => pedigreeMembers[sampleId].SegmentSets[setPosition].SelectedSet = segmentsSet);
            offspringIDs.ForEach(sampleId => pedigreeMembers[sampleId].SegmentSets[setPosition].SelectedSet = segmentsSet);
            var nSegments = pedigreeMembers[parentIDs.First()].SegmentSets[setPosition].GetSet(segmentsSet).Count;
            for (var segmentPosition = 0; segmentPosition < nSegments; segmentPosition++)
            {
                var canvasSegmentIndex = new CanvasSegmentIndex(setPosition, segmentsSet, segmentPosition);

                var copyNumbersLikelihood = AssignCopyNumberWithPedigreeInfo(pedigreeMembers, pedigreeMembersInfo, model, parentIDs,
                    offspringIDs, canvasSegmentIndex, transitionMatrix, offspringsGenotypes);

                EstimateQScoresWithPedigreeInfo(pedigreeMembers, pedigreeMembersInfo, model, parentIDs,
                    offspringIDs, canvasSegmentIndex, copyNumbersLikelihood);

                if (!UseMafInformation(pedigreeMembers, setPosition, segmentPosition, segmentsSet))
                    AssignMccWithPedigreeInfo(pedigreeMembers, pedigreeMembersInfo, model, parentIDs,
                        offspringIDs, canvasSegmentIndex, genotypes);
            }
        }

        private double GetSegmentSetLikelihood(SampleList<PedigreeMember> pedigreeMembers,
            SampleList<PedigreeMemberInfo> pedigreeMembersInfo, SampleList<CopyNumberModel> model, List<SampleId> parentIDs,
            List<SampleId> offspringIDs, int setPosition, SegmentsSet segmentsSet, double[][] transitionMatrix,
            List<List<Genotype>> offspringsGenotypes)
        {
            double segmentSetLikelihood = 0;
            var nSegments = pedigreeMembers[parentIDs.First()].SegmentSets[setPosition].GetSet(segmentsSet).Count;
            for (var segmentPosition = 0; segmentPosition < nSegments; segmentPosition++)
            {
                var index = new CanvasSegmentIndex(setPosition, segmentsSet, segmentPosition);
                segmentSetLikelihood += AssignCopyNumberWithPedigreeInfo(pedigreeMembers, pedigreeMembersInfo, model, parentIDs,
                    offspringIDs, index, transitionMatrix, offspringsGenotypes).MaximalLikelihood;
            }

            segmentSetLikelihood /= nSegments;
            return segmentSetLikelihood;
        }


        /// <summary>
        /// Calculates maximal likelihood for copy numbers. Updated CanvasSegment CopyNumber only. 
        /// </summary>
        /// <param name="samples"></param>
        /// <param name="samplesInfo"></param>
        /// <param name="model"></param>
        /// <param name="parentIDs"></param>
        /// <param name="offspringIDs"></param>
        /// <param name="setPosition"></param>
        /// <param name="segmentPosition"></param>
        /// <param name="segmentsSet"></param>
        /// <param name="transitionMatrix"></param>
        /// <param name="offspringsGenotypes"></param>
        public CopyNumbersLikelihood AssignCopyNumberWithPedigreeInfo(SampleList<PedigreeMember> samples,
            SampleList<PedigreeMemberInfo> samplesInfo, SampleList<CopyNumberModel> model, List<SampleId> parentIDs,
            List<SampleId> offspringIDs, CanvasSegmentIndex index, double[][] transitionMatrix,
                List<List<Genotype>> offspringsGenotypes)
        {
            int nCopies = CallerParameters.MaximumCopyNumber;
            var names = offspringIDs.Union(parentIDs).Select(x => x.ToString()).ToList();
            var density = new CopyNumbersLikelihood(nCopies, names);
            InitializeCn(index, samples);
            density.MaximalLikelihood = 0;
            var coverages = parentIDs.Select(id => Math.Min(samples[id].GetCoverage(index, CallerParameters.NumberOfTrimmedBins),
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
                            double coverage = Math.Min(samples[child].GetCoverage(index, CallerParameters.NumberOfTrimmedBins), samplesInfo[child].MeanCoverage * 3.0);

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
                            samples[parentIDs.First()].GetCanvasSegment(index).CopyNumber = cn1;
                            samples[parentIDs.Last()].GetCanvasSegment(index).CopyNumber = cn2;
                            for (int counter = 0; counter < offspringIDs.Count; counter++)
                            {
                                samples[offspringIDs[counter]].GetCanvasSegment(index).CopyNumber =
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
        public void AssignMccWithPedigreeInfo(SampleList<PedigreeMember> samples,
            SampleList<PedigreeMemberInfo> samplesInfo, SampleList<CopyNumberModel> model, List<SampleId> parentIDs,
            List<SampleId> offspringIDs, CanvasSegmentIndex index, Dictionary<int, List<Genotype>> genotypes)
        {
            double maximalLikelihood = Double.MinValue;
            int parent1CopyNumber = samples[parentIDs.First()].GetCanvasSegment(index).CopyNumber;
            int parent2CopyNumber = samples[parentIDs.Last()].GetCanvasSegment(index).CopyNumber;

            foreach (var parent1GtStates in genotypes[parent1CopyNumber])
            {
                foreach (var parent2GtStates in genotypes[parent2CopyNumber])
                {
                    var bestChildGtStates = new List<Genotype>();
                    double currentLikelihood = 1;
                    foreach (SampleId child in offspringIDs)
                    {
                        int childCopyNumber = samples[child].GetCanvasSegment(index).CopyNumber;
                        bool isInheritedCnv = !samples[child].GetCanvasSegment(index).DqScore.HasValue;
                        double bestLikelihood = Double.MinValue;
                        Genotype bestGtState = null;
                        bestLikelihood = GetProbandLikelihood(model[child], index, genotypes, childCopyNumber,
                            parent1GtStates, parent2GtStates, isInheritedCnv, samples[child], bestLikelihood, ref bestGtState);
                        bestChildGtStates.Add(bestGtState);
                        currentLikelihood *= bestLikelihood;
                    }
                    currentLikelihood *= GetCurrentGtLikelihood(model[parentIDs.First()], samples[parentIDs.First()], index, parent1GtStates) *
                                         GetCurrentGtLikelihood(model[parentIDs.Last()], samples[parentIDs.Last()], index, parent2GtStates);

                    currentLikelihood = Double.IsNaN(currentLikelihood) || Double.IsInfinity(currentLikelihood)
                        ? 0
                        : currentLikelihood;

                    if (currentLikelihood > maximalLikelihood)
                    {
                        maximalLikelihood = currentLikelihood;
                        AssignMCC(samples[parentIDs.First()], model[parentIDs.First()], index, genotypes, parent1GtStates, parent1CopyNumber);
                        AssignMCC(samples[parentIDs.Last()], model[parentIDs.Last()], index, genotypes, parent2GtStates, parent2CopyNumber);
                        var counter = 0;
                        foreach (SampleId child in offspringIDs)
                        {
                            if (bestChildGtStates[counter] == null) continue;
                            int childCopyNumber = samples[child].GetCanvasSegment(index).
                            CopyNumber;
                            AssignMCC(samples[child], model[child], index, genotypes, bestChildGtStates[counter], childCopyNumber);
                            counter++;
                        }
                    }
                }
            }
        }

        private double GetProbandLikelihood(CopyNumberModel copyNumberModel, CanvasSegmentIndex index, Dictionary<int, List<Genotype>> genotypes,
            int childCopyNumber, Genotype parent1GtStates, Genotype parent2GtStates, bool isInheritedCnv, PedigreeMember child,
            double bestLikelihood, ref Genotype bestGtState)
        {
            foreach (var childGtState in genotypes[childCopyNumber])
            {
                double currentChildLikelihood;
                if (IsGtPedigreeConsistent(parent1GtStates, childGtState) &&
                    IsGtPedigreeConsistent(parent2GtStates, childGtState)
                    && isInheritedCnv)
                    currentChildLikelihood = copyNumberModel.GetCurrentGtLikelihood(child.GetAlleleCounts(index), childGtState);
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

        private static void AssignMCC(PedigreeMember pedigreeMember, CopyNumberModel copyNumberModel, CanvasSegmentIndex index,
            Dictionary<int, List<Genotype>> genotypes, Genotype gtStates, int copyNumber)
        {
            const int diploidCopyNumber = 2;
            const int haploidCopyNumber = 1;
            if (copyNumber > diploidCopyNumber)
            {

                pedigreeMember.GetCanvasSegment(index).MajorChromosomeCount =
                    Math.Max(gtStates.CountsA, gtStates.CountsB);
                int? selectedGtState = genotypes[copyNumber].IndexOf(gtStates);
                pedigreeMember.GetCanvasSegment(index).MajorChromosomeCountScore =
                    copyNumberModel.GetGtLikelihoodScore(pedigreeMember.GetAlleleCounts(index), genotypes[copyNumber], ref selectedGtState);
            }
            else
            {
                pedigreeMember.GetCanvasSegment(index).MajorChromosomeCount = copyNumber == diploidCopyNumber
                    ? haploidCopyNumber : copyNumber;
                pedigreeMember.GetCanvasSegment(index).MajorChromosomeCountScore = null;
            }
        }

        private static double GetCurrentGtLikelihood(CopyNumberModel copyNumberModel, PedigreeMember pedigreeMember, CanvasSegmentIndex index, Genotype gtStates)
        {
            return copyNumberModel.GetCurrentGtLikelihood(pedigreeMember.GetAlleleCounts(index), gtStates);
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
        /// <param name="samples"></param>
        /// <param name="segmentPosition"></param>
        /// <param name="copyNumberCombinations"></param>
        public double[][] AssignCopyNumberNoPedigreeInfo(SampleList<PedigreeMember> samples, SampleList<PedigreeMemberInfo> samplesInfo,
            SampleList<CopyNumberModel> copyNumberModel, CanvasSegmentIndex index, List<List<int>> copyNumberCombinations)
        {
            const int defaultCn = 2;
            const double maxCoverageMultiplier = 3.0;

            double maximalLikelihood = 0;
            foreach (PedigreeMember sample in samples.SampleData)
                sample.GetCanvasSegment(index).CopyNumber = defaultCn;
            int nCopies = CallerParameters.MaximumCopyNumber;
            var names = samples.SampleIds.Select(x => x.ToString()).ToList();
            var totalLikelihoods = new List<double>();
            foreach (var copyNumberCombination in copyNumberCombinations)
            {
                double totalLikelihood = 0;
                foreach (var sampleId in samples.SampleIds)
                {
                    maximalLikelihood = 0;
                    foreach (var copyNumber in copyNumberCombination)
                    {
                        var currentLikelihood =
                            copyNumberModel[sampleId].GetCnLikelihood(
                                Math.Min(samples[sampleId].GetCoverage(index, CallerParameters.NumberOfTrimmedBins),
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

            var density = new double[samples.Count()][];
            // no need to iterate over multiple genotypes for n (samples.Count) = 1
            if (samples.Count() == 1)
            {
                samples.Single().Value.GetCanvasSegment(index).CopyNumber =
                    copyNumberCombinations[totalLikelihoods.IndexOf(totalLikelihoods.Max())].First();
                density[0] = totalLikelihoods.ToArray();
                return density;
            }

            var bestcopyNumberCombination = copyNumberCombinations[totalLikelihoods.IndexOf(totalLikelihoods.Max())];
            int counter = 0;
            foreach (var sampleId in samples.SampleIds)
            {
                maximalLikelihood = 0;
                density[counter] = new double[nCopies];
                counter++;
                foreach (var copyNumber in bestcopyNumberCombination)
                {
                    double currentLikelihood =
                        copyNumberModel[sampleId].GetCnLikelihood(
                            Math.Min(samples[sampleId].GetCoverage(index, CallerParameters.NumberOfTrimmedBins),
                                samplesInfo[sampleId].MeanCoverage * maxCoverageMultiplier))[copyNumber];
                    currentLikelihood = Double.IsNaN(currentLikelihood) || Double.IsInfinity(currentLikelihood)
                        ? 0
                        : currentLikelihood;
                    if (currentLikelihood > maximalLikelihood)
                    {
                        maximalLikelihood = currentLikelihood;
                        samples[sampleId].GetCanvasSegment(index).CopyNumber = copyNumber;
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
        public void AssignMccNoPedigreeInfo(SampleList<PedigreeMember> pedigreeMembers, SampleList<PedigreeMemberInfo> pedigreeMembersInfo,
            SampleList<CopyNumberModel> model, CanvasSegmentIndex index, Dictionary<int, List<Genotype>> genotypes)
        {
            foreach (var sampleId in pedigreeMembers.SampleIds)
            {
                int copyNumber = pedigreeMembers[sampleId].GetCanvasSegment(index).CopyNumber;
                if (copyNumber > 2)
                {
                    pedigreeMembers[sampleId].GetCanvasSegment(index).MajorChromosomeCount = copyNumber == 2 ? 1 : copyNumber;
                    return;
                }
                var genotypeset = genotypes[copyNumber];
                int? selectedGtState = null;
                double gqscore = model[sampleId].GetGtLikelihoodScore(pedigreeMembers[sampleId].GetAlleleCounts(index),
                    genotypeset, ref selectedGtState);
                pedigreeMembers[sampleId].GetCanvasSegment(index).MajorChromosomeCountScore = gqscore;
                if (selectedGtState.HasValue)
                    pedigreeMembers[sampleId].GetCanvasSegment(index).MajorChromosomeCount =
                        Math.Max(genotypeset[selectedGtState.Value].CountsA,
                            genotypeset[selectedGtState.Value].CountsB);
            }
        }


        private static void InitializeCn(CanvasSegmentIndex index, SampleList<PedigreeMember> samples)
        {
            const int defaultCn = 2;
            samples.ForEach(sample => sample.Value.GetCanvasSegment(index)
                .CopyNumber = defaultCn);
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


        public static SampleList<PedigreeMember.Kinship> ReadPedigreeFile(string pedigreeFile)
        {
            var kinships = new SampleList<PedigreeMember.Kinship>();
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
                        kinships.Add(new SampleId(fields[1]), PedigreeMember.Kinship.Parent);
                    else if (proband == "affected")
                        kinships.Add(new SampleId(fields[1]), PedigreeMember.Kinship.Proband);
                    else
                        Console.WriteLine($"Unused pedigree member: {row}");
                }
            }
            return kinships;
        }

        public List<double> GetSingleSampleQualityScores(CopyNumbersLikelihood density, List<int> cnStates,
            List<string> sampleNames)
        {
            var singleSampleQualityScores = new List<double>();
            if (density.Count != cnStates.Count)
                throw new ArgumentException("Size of CopyNumbersLikelihood should be equal to number of CN states");
            for (int index = 0; index < sampleNames.Count; index++)
            {
                string sampleName = sampleNames[index];
                var cnMarginalProbabilities = density.GetMarginalProbability(cnStates.Count,
                    CallerParameters.MaximumCopyNumber, sampleName);
                double normalizationConstant = cnMarginalProbabilities.Sum();
                var qscore = -10.0 *
                             Math.Log10((normalizationConstant - cnMarginalProbabilities[cnStates[index]]) /
                                        normalizationConstant);
                if (Double.IsInfinity(qscore) | qscore > CallerParameters.MaxQscore)
                    qscore = CallerParameters.MaxQscore;
                singleSampleQualityScores.Add(qscore);
            }
            return singleSampleQualityScores;
        }

        public double GetConditionalDeNovoQualityScore(CopyNumbersLikelihood density, int probandIndex,
            int probandCopyNumber, string probandName, int parent1Index, int parent2Index,
            List<int> remainingProbandIndex)
        {

            var numerator = 0.0;
            var denominator = 0.0;
            const int diploidState = 2;
            int nSamples = density.Count;
            var probandMarginalProbabilities = density.GetMarginalProbability(nSamples,
                CallerParameters.MaximumCopyNumber, probandName);
            double normalization = probandMarginalProbabilities[probandCopyNumber] +
                                   probandMarginalProbabilities[diploidState];
            double probandMarginalAlt = probandMarginalProbabilities[probandCopyNumber] / normalization;

            foreach (var copyNumberIndex in density.Indices.Where(x => x[probandIndex] == probandCopyNumber))
            {
                if (!(density.GetJointProbability(copyNumberIndex.ToArray()) > 0.0)) continue;

                var holder = density.GetJointProbability(copyNumberIndex);
                denominator += holder;

                if (copyNumberIndex[parent1Index] == diploidState && copyNumberIndex[parent2Index] == diploidState &&
                    remainingProbandIndex.All(index => copyNumberIndex[index] == 2))
                    numerator += holder;
            }

            const double q60 = 0.000001;
            double denovoProbability = (1 - numerator / denominator) * (1 - probandMarginalAlt);
            double qscore = -10.0 * Math.Log10(Math.Max(denovoProbability, q60));
            return qscore;
        }
    }
}

