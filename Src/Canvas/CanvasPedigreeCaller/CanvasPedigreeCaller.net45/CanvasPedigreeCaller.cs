using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection.Metadata;
using System.Runtime.InteropServices;
using MathNet.Numerics.Distributions;
using System.Threading.Tasks;
using CanvasCommon;
using Combinatorics.Collections;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Isas.SequencingFiles;
using MathNet.Numerics;



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

        #endregion

        internal int CallVariantsInPedigree(List<string> variantFrequencyFiles, List<string> segmentFiles,
            string outVcfFile, string ploidyBedPath,
            string referenceFolder, List<string> sampleNames, string commonCNVsbedPath, string pedigreeFile)
        {
            // load files
            // initialize data structures and classes
            int fileCounter = 0;
            var kinships = ReadPedigreeFile(pedigreeFile);
            var pedigreeMembers = new LinkedList<PedigreeMember>();
            foreach (string sampleName in sampleNames)
            {
                var pedigreeMember = SetPedigreeMember(variantFrequencyFiles, segmentFiles, ploidyBedPath, sampleName,
                    fileCounter,
                    CallerParameters.DefaultReadCountsThreshold, referenceFolder, CallerParameters.NumberOfTrimmedBins,
                    commonCNVsbedPath);
                Console.WriteLine($"pedigreeMember for {sampleName} is set");

                pedigreeMember.Kin = kinships[pedigreeMember.Name] == PedigreeMember.Kinship.Parent
                    ? PedigreeMember.Kinship.Parent
                    : PedigreeMember.Kinship.Offspring;
                if (kinships[pedigreeMember.Name] == PedigreeMember.Kinship.Proband)
                {
                    pedigreeMember.Kin = PedigreeMember.Kinship.Proband;
                    pedigreeMembers.AddFirst(pedigreeMember);
                }
                else
                {
                    pedigreeMembers.AddLast(pedigreeMember);

                }

                fileCounter++;
            }

            var maxCoreNumber = 30;
            var segmentIntervals = GetParallelIntervals(pedigreeMembers.First().SegmentSets.Count,
                Math.Min(Environment.ProcessorCount, maxCoreNumber));

            var parents = GetParents(pedigreeMembers);
            var offsprings = GetChildren(pedigreeMembers);
            double[][] transitionMatrix = GetTransitionMatrix(CallerParameters.MaximumCopyNumber);
            var parentalGenotypes = GenerateParentalGenotypes(CallerParameters.MaximumCopyNumber);
            var offspringsGenotypes =
                new List<List<Genotype>>(Convert.ToInt32(Math.Pow(parentalGenotypes.Count, offsprings.Count)));
            GenerateOffspringGenotypes(offspringsGenotypes, parentalGenotypes, offsprings.Count, new List<Genotype>());
            var genotypes = GenerateGenotypeCombinations(CallerParameters.MaximumCopyNumber);

            if (offspringsGenotypes.Count > CallerParameters.MaxNumOffspringGenotypes)
            {
                offspringsGenotypes.Shuffle();
                offspringsGenotypes = offspringsGenotypes.Take(CallerParameters.MaxNumOffspringGenotypes).ToList();
            }

            foreach (PedigreeMember parent in parents)
                parent.CnModel = new CopyNumberModel(CallerParameters.MaximumCopyNumber, parent.MeanCoverage / 2.0,
                    parent.MeanMafCoverage / 2.0, parent.MeanCoverage * 2.5, parent.MeanMafCoverage * 2.5,
                    parent.MaxCoverage);
            foreach (PedigreeMember offspring in offsprings)
                offspring.CnModel = new CopyNumberModel(CallerParameters.MaximumCopyNumber, offspring.MeanCoverage / 2.0,
                    offspring.MeanMafCoverage / 2.0, offspring.MeanCoverage * 2.5, offspring.MeanMafCoverage * 2.5,
                    offspring.MaxCoverage);

            Parallel.ForEach(
                segmentIntervals,
                interval =>
                {
                    Console.WriteLine($"{DateTime.Now} Launching SPW task for segment {interval.Start} - {interval.End}");
                    for (int segmentIndex = interval.Start; segmentIndex <= interval.End; segmentIndex++)
                    {
                        CallVariantInPedigree(pedigreeMembers, parents, offsprings, segmentIndex,
                            transitionMatrix, offspringsGenotypes, genotypes);
                    }
                    Console.WriteLine($"{DateTime.Now} Finished SPW task for segment {interval.Start} - {interval.End}");
                });

            var names = MergeAndExtractSegments(outVcfFile, referenceFolder, pedigreeMembers,
                out List<double?> diploidCoverage,
                out List<List<CanvasSegment>> segments);
            var ploidies = pedigreeMembers.Select(x => x.Ploidy).ToList();
            CanvasSegmentWriter.WriteMultiSampleSegments(outVcfFile, segments, diploidCoverage, referenceFolder, names,
                null, ploidies,
                QualityFilterThreshold, isPedigreeInfoSupplied: true,
                denovoQualityThreshold: DeNovoQualityFilterThreshold);

            var outputFolder = new FileLocation(outVcfFile).Directory;
            foreach (var pedigreeMember in pedigreeMembers)
            {
                var outputVcfPath = SingleSampleCallset.GetSingleSamplePedigreeVcfOutput(outputFolder,
                    pedigreeMember.Name);
                CanvasSegmentWriter.WriteSegments(outputVcfPath.FullName, pedigreeMember.Segments,
                    pedigreeMember.MeanCoverage, referenceFolder,
                    pedigreeMember.Name, null, pedigreeMember.Ploidy, QualityFilterThreshold,
                    isPedigreeInfoSupplied: true,
                    denovoQualityThreshold: DeNovoQualityFilterThreshold);
            }
            return 0;
        }

        /// <summary>
        /// Derives metrics from b-allele counts within each segment and determines whereas to use them for calculating MCC
        /// </summary>
        /// <param name="pedigreeMembers"></param>
        /// <param name="segmentIndex"></param>
        /// <returns></returns>
        private bool UseMafInformation(LinkedList<PedigreeMember> pedigreeMembers, int segmentSetIndex, int segmentIndex, SegmentsSet segmentsSet)
        {
            const double copyNumberCoverageOneAndHalfNormalizer = 0.75;
            var alleles = pedigreeMembers.Select(
                x => x.SegmentSets[segmentSetIndex].GetSet(segmentsSet)[segmentIndex].Balleles?
                    .BAlleles.Select(y => y.Counts).ToList());
            var alleleCounts = alleles.Select(allele => allele?.Count ?? 0).ToList();
            bool lowAlleleCounts = alleleCounts.Select(x => x < CallerParameters.DefaultReadCountsThreshold).Any(c => c == true);
            var coverageCounts = pedigreeMembers.Select(x => x.SegmentSets[segmentSetIndex].GetSet(segmentsSet)[segmentIndex].MedianCount).ToList();
            var isSkewedHetHomRatio = false;
            if (false)
            {
                var minorMeanMAF =
                    alleles.Select(allele => allele.Select(y => Math.Min(y.Item1, y.Item2)).ToList().Average());
                var minorMedianMAF = pedigreeMembers.Select(x =>x.SegmentSets[segmentSetIndex].GetSet(segmentsSet)[segmentIndex]?.Balleles.MedianCounts.Item2);
                var isHighCoverageCounts = pedigreeMembers.Select(x => x.MeanCoverage).Zip(coverageCounts,
                    (meanCoverage, segmenCoverage) =>
                        segmenCoverage > meanCoverage * copyNumberCoverageOneAndHalfNormalizer);
                var zip = minorMeanMAF.Zip(minorMedianMAF, (mn, md) => new {mn, md})
                    .Zip(isHighCoverageCounts, (t, cov) =>
                        new {meanMAF = t.mn, medianMAF = t.md, isHighCoverage = cov});
                isSkewedHetHomRatio = zip.Any(sample => sample.meanMAF > sample.medianMAF && sample.medianMAF == 0 &&
                                                        sample.isHighCoverage);
            }
            double alleleDensity = pedigreeMembers.First().SegmentSets[segmentSetIndex].GetSet(segmentsSet)[segmentIndex].Length /
                                   Math.Max(alleleCounts.Average(), 1.0);
            bool useCnLikelihood = lowAlleleCounts ||
                                   alleleDensity < CallerParameters.DefaultAlleleDensityThreshold ||
                                   alleleCounts.Any(x => x > CallerParameters.DefaultPerSegmentAlleleMaxCounts) ||
                                   coverageCounts.Any(coverage => coverage < MedianCoverageThreshold) ||
                                   isSkewedHetHomRatio;
            // for now only use lowAlleleCounts metric
            return lowAlleleCounts;
        }

        /// <summary>
        /// Merges segments and extracts coverage / segments and pedigreeMembers
        /// </summary>
        /// <param name="outVcfFile"></param>
        /// <param name="referenceFolder"></param>
        /// <param name="pedigreeMembers"></param>
        /// <param name="meanCoverage"></param>
        /// <param name="segments"></param>
        /// <returns></returns>
        private List<string> MergeAndExtractSegments(string outVcfFile, string referenceFolder,
            LinkedList<PedigreeMember> pedigreeMembers, out List<double?> meanCoverage,
            out List<List<CanvasSegment>> segments)
        {
            MergeSegments(pedigreeMembers, CallerParameters.MinimumCallSize);
            var names = pedigreeMembers.Select(x => x.Name).ToList();
            meanCoverage = pedigreeMembers.Select(x => (double?) x.MeanCoverage).ToList();
            segments = pedigreeMembers.Select(x => x.Segments).ToList();

            var outputFolder = new FileLocation(outVcfFile).Directory;
            foreach (var member in pedigreeMembers)
            {
                var coverageOutputPath = SingleSampleCallset.GetCoverageAndVariantFrequencyOutput(outputFolder,
                    member.Name);
                CanvasSegment.WriteCoveragePlotData(member.Segments, member.MeanCoverage, member.Ploidy, coverageOutputPath.FullName, referenceFolder);
            }
            return names;
        }

        internal int CallVariants(List<string> variantFrequencyFiles, List<string> segmentFiles, string outVcfFile,
            string ploidyBedPath, string referenceFolder, List<string> sampleNames, string commonCNVsbedPath)
        {
            // load files
            // initialize data structures and classes
            var fileCounter = 0;
            var pedigreeMembers = new LinkedList<PedigreeMember>();
            foreach (string sampleName in sampleNames)
            {
                var pedigreeMember = SetPedigreeMember(variantFrequencyFiles, segmentFiles, ploidyBedPath, sampleName,
                    fileCounter,
                    CallerParameters.DefaultReadCountsThreshold, referenceFolder, CallerParameters.NumberOfTrimmedBins,
                    commonCNVsbedPath);
                pedigreeMembers.AddLast(pedigreeMember);
                fileCounter++;
            }

            const int maxCoreNumber = 30;
            var segmentIntervals = GetParallelIntervals(pedigreeMembers.First().SegmentSets.Count,
                Math.Min(Environment.ProcessorCount, maxCoreNumber));
            var genotypes = GenerateGenotypeCombinations(CallerParameters.MaximumCopyNumber);
            int maxAlleleNumber = Math.Min(CallerParameters.MaxAlleleNumber, pedigreeMembers.Count);
            var copyNumberCombinations = GenerateCopyNumberCombinations(CallerParameters.MaximumCopyNumber,
                maxAlleleNumber);

            foreach (PedigreeMember pedigreeMember in pedigreeMembers)
                pedigreeMember.CnModel = new CopyNumberModel(CallerParameters.MaximumCopyNumber,
                    pedigreeMember.MeanCoverage / 2.0, pedigreeMember.MeanMafCoverage / 2.0,
                    pedigreeMember.Variance, pedigreeMember.MafVariance, pedigreeMember.MaxCoverage);

            Parallel.ForEach(
                segmentIntervals,
                interval =>
                {
                    Console.WriteLine($"{DateTime.Now} Launching SPW task for segment {interval.Start} - {interval.End}");
                    for (int segmentIndex = interval.Start; segmentIndex <= interval.End; segmentIndex++)
                    {
                        CallVariant(pedigreeMembers, segmentIndex, copyNumberCombinations, genotypes);
                    }
                    Console.WriteLine($"{DateTime.Now} Finished SPW task for segment {interval.Start} - {interval.End}");
                });

            var names = MergeAndExtractSegments(outVcfFile, referenceFolder, pedigreeMembers,
                out List<double?> diploidCoverage,
                out List<List<CanvasSegment>> segments);
            var ploidies = pedigreeMembers.Select(x => x.Ploidy).ToList();

            CanvasSegmentWriter.WriteMultiSampleSegments(outVcfFile, segments, diploidCoverage, referenceFolder, names,
                null, ploidies,
                QualityFilterThreshold, isPedigreeInfoSupplied: false);

            var outputFolder = new FileLocation(outVcfFile).Directory;
            foreach (var pedigreeMember in pedigreeMembers)
            {
                var outputVcfPath = SingleSampleCallset.GetSingleSamplePedigreeVcfOutput(outputFolder,
                    pedigreeMember.Name);
                CanvasSegmentWriter.WriteSegments(outputVcfPath.FullName, pedigreeMember.Segments,
                    pedigreeMember.MeanCoverage, referenceFolder,
                    pedigreeMember.Name, null, pedigreeMember.Ploidy, QualityFilterThreshold,
                    isPedigreeInfoSupplied: false);
            }
            return 0;
        }


        private static void MergeSegments(LinkedList<PedigreeMember> pedigreeMembers, int minimumCallSize)
        {
            foreach (var pedigreeMember in pedigreeMembers) { 
                pedigreeMember.Segments = pedigreeMember.SegmentSets.SelectMany(x => x.GetSet(x.SelectedSet)).ToList();
                pedigreeMember.Segments = pedigreeMember.Segments.OrderBy(c1 => c1.Chr).ThenBy(c2 => c2.Begin).ToList();
            }

            int nSegments = pedigreeMembers.First().Segments.Count;
            var copyNumbers = new List<List<int>>(nSegments);
            var qscores = new List<double>(nSegments);
            foreach (int segmentIndex in Enumerable.Range(0, nSegments))
            {
                copyNumbers.Add(new List<int>());
                qscores.Add(0);
                foreach (var pedigreeMember in pedigreeMembers)
                {
                    copyNumbers[segmentIndex].Add(pedigreeMember.Segments[segmentIndex].CopyNumber);
                    qscores[segmentIndex] += pedigreeMember.Segments[segmentIndex].QScore;
                }
                qscores[segmentIndex] /= pedigreeMembers.Count;
            }

            if (copyNumbers == null && qscores != null || copyNumbers != null & qscores == null)
                throw new ArgumentException("Both copyNumbers and qscores arguments must be specified.");
            if (copyNumbers != null && copyNumbers.Count != pedigreeMembers.First().Segments.Count)
                throw new ArgumentException("Length of copyNumbers list should be equal to the number of segments.");
            if (qscores != null && qscores.Count != pedigreeMembers.First().Segments.Count)
                throw new ArgumentException("Length of qscores list should be equal to the number of segments.");

            foreach (var pedigreeMember in pedigreeMembers)
                CanvasSegment.MergeSegments(ref pedigreeMember.Segments, minimumCallSize, 10000, copyNumbers, qscores);
        }

        private static PedigreeMember SetPedigreeMember(List<string> variantFrequencyFiles, List<string> segmentFiles,
            string ploidyBedPath,
            string sampleName, int fileCounter, int defaultAlleleCountThreshold, string referenceFolder,
            int numberOfTrimmedBins,
            string commonCNVsbedPath)
        {
            var pedigreeMember = new PedigreeMember {Name = sampleName};
            var segments = CanvasSegment.ReadSegments(segmentFiles[fileCounter]);
            pedigreeMember.MeanMafCoverage = CanvasIO.LoadFrequenciesBySegment(variantFrequencyFiles[fileCounter],
                segments, referenceFolder);
            foreach (var segment in segments)
                if (segment.Balleles.BAlleles.Count > defaultAlleleCountThreshold)
                    segment.Balleles.MedianCounts = Balleles.SetMedianCounts(segment.Balleles);
            pedigreeMember.Variance = GetCoverageVariance(numberOfTrimmedBins, segments);
            pedigreeMember.MafVariance = GetMafVariance(segments);
            pedigreeMember.MeanCoverage = segments.Any()
                ? segments.Select(x => x.MedianCount).Average()
                : 0;
            pedigreeMember.MaxCoverage = segments.Any()
                ? (int) (segments.Select(x => x.MedianCount).Max() + 10)
                : 0;
            if (!ploidyBedPath.IsNullOrEmpty() && File.Exists(ploidyBedPath))
                pedigreeMember.Ploidy = PloidyInfo.LoadPloidyFromVcfFile(ploidyBedPath, pedigreeMember.Name);

            if (commonCNVsbedPath != null)
            {
                var segmentsByChromosome = CanvasSegment.GetSegmentsByChromosome(segments);
                var coverage = CanvasSegment.ReadBEDInput(segmentFiles[fileCounter]);
                var commonRegions = CanvasCommon.Utilities.LoadBedFile(commonCNVsbedPath);
                CanvasCommon.Utilities.SortAndOverlapCheck(commonRegions, commonCNVsbedPath);
                if (IdenticalChromosomeNames(commonRegions, coverage) == 0)
                    throw new ArgumentException(
                        $"Chromosome names in a common CNVs bed file {commonCNVsbedPath} does not match " +
                        $"chromosomes in {segmentFiles[fileCounter]}");

                var segmentIntervalsByChromosome = new Dictionary<string, List<Interval>>();
                Parallel.ForEach(commonRegions.Keys, chr => segmentIntervalsByChromosome[chr] = 
                CanvasSegment.RemapCommonRegions(commonRegions[chr],coverage.StartByChr[chr], coverage.EndByChr[chr]));
                var allelesByChromosome = CanvasIO.ReadFrequencies(variantFrequencyFiles[fileCounter], segmentIntervalsByChromosome,
                    referenceFolder, out float meanCoverage);

                var segmentsSetByChromosome = new Dictionary<string, List<CanvasSegmentsSet>>();
                Parallel.ForEach(
                    segmentsByChromosome.Keys,
                    chr =>
                    {
                        if (commonRegions.Keys.Any(chromosome => chromosome == chr)) {
                            Console.WriteLine($"CreateSegmentsFromCommonCnvs for {chr} ");
                            var commonCnvCanvasSegments = CanvasSegment.CreateSegmentsFromCommonCnvs(coverage, chr,
                                segmentIntervalsByChromosome[chr]);
                            for (var index = 0; index < commonCnvCanvasSegments.Count; index++)
                            {
                                foreach (var genotype in allelesByChromosome[chr][index])
                                    commonCnvCanvasSegments[index].Balleles.BAlleles.Add(CanvasIO.GetAllele(genotype));
                                if (commonCnvCanvasSegments[index].Balleles.BAlleles.Count > defaultAlleleCountThreshold)
                                    commonCnvCanvasSegments[index].Balleles.MedianCounts = Balleles.SetMedianCounts(commonCnvCanvasSegments[index].Balleles);
                            }
                            segmentsSetByChromosome[chr] = CanvasSegment.MergeCommonCnvSegments(segmentsByChromosome[chr], commonCnvCanvasSegments, chr, defaultAlleleCountThreshold) ?? 
                            segmentsByChromosome[chr].Select(segment => new CanvasSegmentsSet(setA: new List<CanvasSegment> {segment}, setB: null)).ToList();
                            Console.WriteLine($"SegmentsFromCommonCnvs for {chr} returned");
                        }
                        else
                        {
                            segmentsSetByChromosome[chr] = segmentsByChromosome[chr].Select(segment => 
                            new CanvasSegmentsSet(setA: new List<CanvasSegment> { segment }, setB: null)).ToList();
                        }
                    });
                Console.WriteLine($"Create SegmentSets for {sampleName} ");
                pedigreeMember.SegmentSets.AddRange(segmentsSetByChromosome.OrderBy(i => i.Key).Select(x => x.Value).SelectMany(x=>x).ToList());
            }
            else
            {
                Console.WriteLine($"Alternative");

                pedigreeMember.SegmentSets =
                    segments.Select(
                            segment =>
                                new CanvasSegmentsSet(setA: new List<CanvasSegment> {segment}, setB: null))
                        .ToList();
            }
            Console.WriteLine($"Done");
            return pedigreeMember;
        }


        private static double GetCoverageVariance(int numberOfTrimmedBins, List<CanvasSegment> segments)
        {
            return
                Math.Pow(
                    Utilities.StandardDeviation(
                        segments.Select(x => x.MedianCount).ToArray()), 2);
        }

        private static double GetMafVariance(List<CanvasSegment> segments)
        {
            return Math.Pow(
                Utilities.StandardDeviation(
                    segments.Where(x => x.Balleles.TotalCoverage.Count > 0)
                        .Select(x => x.Balleles.TotalCoverage.Average())
                        .ToArray()), 2);
        }

        private static int IdenticalChromosomeNames(Dictionary<string, List<SampleGenomicBin>> commonRegions,
            CoverageInfo coverage)
        {
            var chromsomes = new HashSet<string>(coverage.CoverageByChr.Keys);
            return commonRegions.Keys.Count(chromosome => chromsomes.Contains(chromosome));
        }

        private void EstimateQScoresWithPedigreeInfo(List<PedigreeMember> parents, List<PedigreeMember> offsprings,
            int segmentSetIndex,
            int segmentIndex, SegmentsSet segmentsSet, CopyNumberDistribution copyNumberLikelihoods)
        {
            var cnStates = GetCnStates(parents, offsprings, segmentSetIndex, segmentIndex, segmentsSet,
                CallerParameters.MaximumCopyNumber);
            var names = parents.Select(x => x.Name).Concat(offsprings.Select(x => x.Name)).ToList();
            var probands = GetProbands(offsprings);
            var singleSampleQualityScores = GetSingleSampleQualityScores(copyNumberLikelihoods, cnStates, names);
            var counter = 0;
            foreach (var sample in parents.Concat(offsprings))
            {
                sample.SegmentSets[segmentSetIndex].GetSet(segmentsSet)[segmentIndex].QScore =
                    singleSampleQualityScores[counter];
                if (sample.SegmentSets[segmentSetIndex].GetSet(segmentsSet)[segmentIndex].QScore <
                    QualityFilterThreshold)
                    sample.SegmentSets[segmentSetIndex].GetSet(segmentsSet)[segmentIndex].Filter =
                        $"q{QualityFilterThreshold}";
                counter++;
            }
            SetDenovoQualityScores(parents, segmentSetIndex, segmentIndex, segmentsSet, copyNumberLikelihoods,
                names, probands, cnStates, singleSampleQualityScores);
        }

        private void SetDenovoQualityScores(List<PedigreeMember> parents, int segmentSetIndex, int segmentIndex,
            SegmentsSet segmentsSet,
            CopyNumberDistribution copyNumberLikelihoods, List<string> names, List<PedigreeMember> probands,
            List<int> cnStates, List<double> singleSampleQualityScores)
        {
            int parent1Index = names.IndexOf(parents.First().Name);
            int parent2Index = names.IndexOf(parents.Last().Name);

            foreach (var proband in probands)
            {
                int probandIndex = names.IndexOf(proband.Name);
                var remainingProbandIndex = probands.Except(proband.ToEnumerable()).Select(x => names.IndexOf(x.Name));
                if (cnStates[probandIndex] != proband.GetPloidy(segmentSetIndex, segmentIndex, segmentsSet) &&
                    // targeted proband is ALT
                    (ParentsRefCheck(parents, segmentSetIndex, segmentIndex, segmentsSet, cnStates, parent1Index,
                         parent2Index) ||
                     // either parent are REF or 
                     IsNotCommonCnv(parents, proband, cnStates, parent1Index, parent2Index, probandIndex, segmentSetIndex,
                         segmentIndex, segmentsSet)) &&
                    // or a common variant 
                    remainingProbandIndex.All(
                        index =>
                            cnStates[index] ==
                            probands[index].GetPloidy(segmentSetIndex, segmentIndex, segmentsSet) ||
                            IsNotCommonCnv(parents, probands[index], cnStates, parent1Index,
                                parent2Index, index, segmentSetIndex, segmentIndex, segmentsSet)) &&
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
                    proband.SegmentSets[segmentSetIndex].GetSet(segmentsSet)[segmentIndex].DqScore =
                        deNovoQualityScore;
                }
            }
        }

        private static bool IsNotCommonCnv(List<PedigreeMember> parents, PedigreeMember proband, List<int> cnStates,
            int parent1Index, int parent2Index,
            int probandIndex, int segmentSetIndex, int segmentIndex, SegmentsSet segmentsSet)
        {
            var parent1Genotypes = GenerateCnAlleles(cnStates[parent1Index]);
            var parent2Genotypes = GenerateCnAlleles(cnStates[parent2Index]);
            var probandGenotypes = GenerateCnAlleles(cnStates[probandIndex]);

            bool isCommoCnv = (parent1Genotypes.Intersect(probandGenotypes).Any() &&
                               parents.First().GetPloidy(segmentSetIndex, segmentIndex, segmentsSet) ==
                               proband.GetPloidy(segmentSetIndex, segmentIndex, segmentsSet)) ||
                              (parent2Genotypes.Intersect(probandGenotypes).Any() &&
                               parents.Last().GetPloidy(segmentSetIndex, segmentIndex, segmentsSet) ==
                               proband.GetPloidy(segmentSetIndex, segmentIndex, segmentsSet));
            return !isCommoCnv;
        }

        private static bool ParentsRefCheck(List<PedigreeMember> parents, int segmentSetIndex, int segmentIndex,
            SegmentsSet segmentsSet,
            List<int> cnStates, int parent1Index, int parent2Index)
        {
            return cnStates[parent1Index] == parents.First().GetPloidy(segmentSetIndex, segmentIndex, segmentsSet) &&
                   cnStates[parent2Index] == parents.Last().GetPloidy(segmentSetIndex, segmentIndex, segmentsSet);
        }

        private void EstimateQScoresNoPedigreeInfo(LinkedList<PedigreeMember> samples, int setPosition,
                int segmentPosition, SegmentsSet segmentsSet, double[][] copyNumberLikelihoods)
        {
            var cnStates =
                samples.Select(
                    x => Math.Min(x.SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].CopyNumber, CallerParameters.MaximumCopyNumber - 1)).ToList();
            int counter = 0;
            foreach (PedigreeMember sample in samples)
            {
                double normalizationConstant = copyNumberLikelihoods[counter].Sum();
                double qscore = -10.0 *
                                Math.Log10((normalizationConstant - copyNumberLikelihoods[counter][cnStates[counter]]) /
                                           normalizationConstant);
                if (Double.IsInfinity(qscore) | qscore > CallerParameters.MaxQscore)
                    qscore = CallerParameters.MaxQscore;
                sample.SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].QScore = qscore;
                if (sample.SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].QScore < QualityFilterThreshold)
                    sample.SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].Filter = $"q{QualityFilterThreshold}";
                counter++;
            }
        }


        private static List<int> GetCnStates(IEnumerable<PedigreeMember> parents, IEnumerable<PedigreeMember> offsprings,
            int segmentSetIndex, int segmentIndex, SegmentsSet segmentsSet, int maximumCopyNumber)
        {
            return
                parents.Select(
                        x =>
                            Math.Min(
                                x.SegmentSets[segmentSetIndex].GetSet(segmentsSet)[segmentIndex]
                                    .CopyNumber, maximumCopyNumber - 1))
                    .Concat(offsprings.Select(
                        x =>
                            Math.Min(
                                x.SegmentSets[segmentSetIndex].GetSet(segmentsSet)[segmentIndex]
                                    .CopyNumber, maximumCopyNumber - 1))).ToList();
        }

        private static List<PedigreeMember> GetChildren(IEnumerable<PedigreeMember> pedigreeMembers)
        {
            return pedigreeMembers.Select(x => x).Where(x => x.Kin != PedigreeMember.Kinship.Parent).ToList();
        }

        private static List<PedigreeMember> GetParents(IEnumerable<PedigreeMember> pedigreeMembers)
        {
            return pedigreeMembers.Select(x => x).Where(x => x.Kin == PedigreeMember.Kinship.Parent).ToList();
        }

        private static List<PedigreeMember> GetProbands(IEnumerable<PedigreeMember> pedigreeMembers)
        {
            return pedigreeMembers.Select(x => x).Where(x => x.Kin == PedigreeMember.Kinship.Proband).ToList();
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
        public void CallVariant(LinkedList<PedigreeMember> samples, int setPosition,
            List<List<int>> copyNumbers, Dictionary<int, List<Genotype>> genotypes)
        {
            SegmentsSet segmentsSet;

            if (samples.First().SegmentSets[setPosition].SetA == null)
                segmentsSet = SegmentsSet.SetB;

            else if (samples.First().SegmentSets[setPosition].SetB == null)
                segmentsSet = SegmentsSet.SetA;

            else
            {
                double setAlikelihood = 0;
                if (samples.First().SegmentSets[setPosition].SetA != null)
                {
                    for (var segmentPosition = 0;
                        segmentPosition < samples.First().SegmentSets[setPosition].SetA.Count;
                        segmentPosition++)
                        setAlikelihood += Utilities.MaxValue(
                            AssignCopyNumberNoPedigreeInfo(samples, setPosition, segmentPosition,
                                SegmentsSet.SetA, copyNumbers));
                    setAlikelihood /= samples.First().SegmentSets[setPosition].SetA.Count;

                }

                double setBlikelihood = 0;
                if (samples.First().SegmentSets[setPosition].SetB != null)
                {
                    for (var segmentPosition = 0;
                        segmentPosition < samples.First().SegmentSets[setPosition].SetB.Count;
                        segmentPosition++)
                        setBlikelihood += Utilities.MaxValue(
                            AssignCopyNumberNoPedigreeInfo(samples, setPosition, segmentPosition,
                                SegmentsSet.SetB, copyNumbers));
                    setBlikelihood /= samples.First().SegmentSets[setPosition].SetB.Count;
                }

                segmentsSet = setAlikelihood > setBlikelihood
                    ? SegmentsSet.SetA
                    : SegmentsSet.SetB;
            }

            samples.ForEach(sample => sample.SegmentSets[setPosition].SelectedSet = segmentsSet);

            for (var segmentPosition = 0;
                segmentPosition <
                samples.First().SegmentSets[setPosition].GetSet(segmentsSet).Count;
                segmentPosition++)
            {
                var ll = AssignCopyNumberNoPedigreeInfo(samples, setPosition, segmentPosition, segmentsSet, copyNumbers);
                EstimateQScoresNoPedigreeInfo(samples, setPosition, segmentPosition, segmentsSet, ll);
                AssignMccNoPedigreeInfo(samples, setPosition, segmentPosition, segmentsSet, genotypes);
            }
        }

        /// <summary>
        /// Identify variant with the highest likelihood at a given setPosition and assign relevant scores
        /// </summary>
        /// <param name="pedigreeMembers"></param>
        /// <param name="parents"></param>
        /// <param name="children"></param>
        /// <param name="setPosition"></param>
        /// <param name="transitionMatrix"></param>
        /// <param name="offspringsGenotypes"></param>
        /// <param name="genotypes"></param>
        public void CallVariantInPedigree(LinkedList<PedigreeMember> pedigreeMembers,
                List<PedigreeMember> parents, List<PedigreeMember> children, int setPosition, double[][] transitionMatrix, List<List<Genotype>> offspringsGenotypes,
                Dictionary<int, List<Genotype>> genotypes)
        {
            SegmentsSet segmentsSet;

            if (parents.First().SegmentSets[setPosition].SetA == null)
                segmentsSet = SegmentsSet.SetB;

            else if (parents.First().SegmentSets[setPosition].SetB == null)
                segmentsSet = SegmentsSet.SetA;

            else
            {
                double setAlikelihood = 0;
                if (parents.First().SegmentSets[setPosition].SetA != null)
                {
                    for (var segmentPosition = 0;
                        segmentPosition < parents.First().SegmentSets[setPosition].SetA.Count;
                        segmentPosition++)
                        setAlikelihood +=
                            AssignCopyNumberWithPedigreeInfo(parents, children, setPosition, segmentPosition,
                                SegmentsSet.SetA, transitionMatrix, offspringsGenotypes).MaximalLikelihood;
                    setAlikelihood /= parents.First().SegmentSets[setPosition].SetA.Count;

                }


                double setBlikelihood = 0;
                if (parents.First().SegmentSets[setPosition].SetB != null)
                {
                    for (var segmentPosition = 0;
                        segmentPosition < parents.First().SegmentSets[setPosition].SetB.Count;
                        segmentPosition++)
                        setBlikelihood +=
                            AssignCopyNumberWithPedigreeInfo(parents, children, setPosition, segmentPosition,
                                SegmentsSet.SetB, transitionMatrix, offspringsGenotypes).MaximalLikelihood;
                    setBlikelihood /= parents.First().SegmentSets[setPosition].SetB.Count;
                }

                segmentsSet = setAlikelihood > setBlikelihood
                    ? SegmentsSet.SetA
                    : SegmentsSet.SetB;
            }

            parents.ForEach(sample => sample.SegmentSets[setPosition].SelectedSet = segmentsSet);
            children.ForEach(sample => sample.SegmentSets[setPosition].SelectedSet = segmentsSet);

            for (var segmentPosition = 0;
                segmentPosition <
                parents.First().SegmentSets[setPosition].GetSet(segmentsSet).Count;
                segmentPosition++)
            {
                var ll = AssignCopyNumberWithPedigreeInfo(parents, children, setPosition, segmentPosition,
                    segmentsSet, transitionMatrix, offspringsGenotypes);
                EstimateQScoresWithPedigreeInfo(parents, children, setPosition, segmentPosition,
                    segmentsSet, ll);
                if (!UseMafInformation(pedigreeMembers, setPosition, segmentPosition,
                    segmentsSet))
                    AssignMccWithPedigreeInfo(parents, children, setPosition, segmentPosition, segmentsSet, genotypes);
            }
        }


        /// <summary>
        /// Calculates maximal likelihood for copy numbers. Updated CanvasSegment CopyNumber only. 
        /// </summary>
        /// <param name="parents"></param>
        /// <param name="children"></param>
        /// <param name="segmentPosition"></param>
        /// <param name="transitionMatrix"></param>
        public CopyNumberDistribution AssignCopyNumberWithPedigreeInfo(List<PedigreeMember> parents,
                List<PedigreeMember> children, int setPosition,
                int segmentPosition, SegmentsSet segmentsSet, double[][] transitionMatrix,
                List<List<Genotype>> offspringsGenotypes)
            {
                int nCopies = CallerParameters.MaximumCopyNumber;
                var names = parents.Select(x => x.Name).Union(children.Select(x => x.Name)).ToList();
                var density = new CopyNumberDistribution(nCopies, names);
                InitializeCn(setPosition, segmentPosition, segmentsSet, parents, children);
                density.MaximalLikelihood = 0;
                var parent1Likelihood =
                    parents.First()
                        .CnModel.GetCnLikelihood(
                            Math.Min(
                                parents.First()
                                    .GetCoverage(setPosition, segmentPosition, segmentsSet,
                                        CallerParameters.NumberOfTrimmedBins), parents.First().MeanCoverage * 3.0));
                var parent2Likelihood =
                    parents.Last()
                        .CnModel.GetCnLikelihood(
                            Math.Min(
                                parents.Last()
                                    .GetCoverage(setPosition, segmentPosition, segmentsSet,
                                        CallerParameters.NumberOfTrimmedBins), parents.Last().MeanCoverage * 3.0));

                if (parent1Likelihood.Count != parent2Likelihood.Count)
                    throw new ArgumentException("Both parents should have the same number of CN states");

                for (int cn1 = 0; cn1 < nCopies; cn1++)
                {
                    for (int cn2 = 0; cn2 < nCopies; cn2++)
                    {
                        foreach (var offspringGtStates in offspringsGenotypes)
                        {
                            double currentLikelihood = parent1Likelihood[cn1] * parent2Likelihood[cn2];
                            int counter = 0;
                            foreach (PedigreeMember child in children)
                            {
                                var modelIndex =
                                    Math.Min(offspringGtStates[counter].CountsA + offspringGtStates[counter].CountsB,
                                        CallerParameters.MaximumCopyNumber - 1);
                                currentLikelihood *= transitionMatrix[cn1][offspringGtStates[counter].CountsA] *
                                                     transitionMatrix[cn2][offspringGtStates[counter].CountsB] *
                                                     child.CnModel.GetCnLikelihood(child.GetCoverage(setPosition,
                                                         segmentPosition, segmentsSet,
                                                         CallerParameters.NumberOfTrimmedBins))[modelIndex];
                                counter++;

                            }
                            int[] copyNumberIndices = {cn1, cn2};
                            var index =
                                copyNumberIndices.Concat(offspringGtStates.Select(x => x.CountsA + x.CountsB)).ToArray();
                            density.SetJointProbability(
                                Math.Max(currentLikelihood, density.GetJointProbability(index)), index);

                            currentLikelihood = Double.IsNaN(currentLikelihood) || Double.IsInfinity(currentLikelihood)
                                ? 0
                                : currentLikelihood;

                            if (currentLikelihood > density.MaximalLikelihood)
                            {
                                density.MaximalLikelihood = currentLikelihood;
                                parents.First().SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].CopyNumber = cn1;
                                parents.Last().SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].CopyNumber = cn2;
                                counter = 0;
                                foreach (PedigreeMember child in children)
                                {
                                    child.SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].CopyNumber = offspringGtStates[counter].CountsA +
                                                                                 offspringGtStates[counter].CountsB;
                                    counter++;
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
            /// <param name="parents"></param>
            /// <param name="children"></param>
            /// <param name="segmentPosition"></param>
            /// <param name="genotypes"></param>
            public void AssignMccWithPedigreeInfo(List<PedigreeMember> parents, List<PedigreeMember> children,
                int setPosition, int segmentPosition, SegmentsSet segmentsSet, Dictionary<int, List<Genotype>> genotypes)
            {
                double maximalLikelihood = Double.MinValue;
                int parent1CopyNumber = parents.First().SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].CopyNumber;
                int parent2CopyNumber = parents.Last().SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].CopyNumber;

                foreach (var parent1GtStates in genotypes[parent1CopyNumber])
                {
                    foreach (var parent2GtStates in genotypes[parent2CopyNumber])
                    {
                        var bestChildGtStates = new List<Genotype>();
                        double currentLikelihood = 1;
                        foreach (PedigreeMember child in children)
                        {
                            int childCopyNumber = child.SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].CopyNumber;
                            bool isInheritedCnv = !child.SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].DqScore.HasValue;
                            double bestLikelihood = Double.MinValue;
                            Genotype bestGtState = null;
                            bestLikelihood = GetProbandLikelihood(setPosition, segmentPosition, segmentsSet, genotypes, childCopyNumber,
                                parent1GtStates, parent2GtStates, isInheritedCnv, child, bestLikelihood, ref bestGtState);
                            bestChildGtStates.Add(bestGtState);
                            currentLikelihood *= bestLikelihood;
                        }
                        currentLikelihood *= GetCurrentGtLikelihood(parents.First(), setPosition, segmentPosition, segmentsSet, parent1GtStates) *
                                             GetCurrentGtLikelihood(parents.Last(), setPosition, segmentPosition, segmentsSet, parent2GtStates);

                        currentLikelihood = Double.IsNaN(currentLikelihood) || Double.IsInfinity(currentLikelihood)
                            ? 0
                            : currentLikelihood;

                        if (currentLikelihood > maximalLikelihood)
                        {
                            maximalLikelihood = currentLikelihood;
                            AssignMCC(parents.First(), setPosition, segmentPosition, segmentsSet, genotypes, parent1GtStates, parent1CopyNumber);
                            AssignMCC(parents.Last(), setPosition, segmentPosition, segmentsSet, genotypes, parent2GtStates, parent2CopyNumber);
                            var counter = 0;
                            foreach (PedigreeMember child in children)
                            {
                                if (bestChildGtStates[counter] == null) continue;
                                int childCopyNumber = child.SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].
                                CopyNumber;
                                AssignMCC(child, setPosition, segmentPosition, segmentsSet, genotypes, bestChildGtStates[counter], childCopyNumber);
                                counter++;
                            }
                        }
                    }
                }
            }

            private double GetProbandLikelihood(int setPosition, int segmentPosition, SegmentsSet segmentsSet, Dictionary<int, List<Genotype>> genotypes,
                int childCopyNumber, Genotype parent1GtStates, Genotype parent2GtStates, bool isInheritedCnv, PedigreeMember child,
                double bestLikelihood, ref Genotype bestGtState)
            {
                foreach (var childGtState in genotypes[childCopyNumber])
                {
                    double currentChildLikelihood;
                    if (IsGtPedigreeConsistent(parent1GtStates, childGtState) &&
                        IsGtPedigreeConsistent(parent2GtStates, childGtState)
                        && isInheritedCnv)
                        currentChildLikelihood = child.CnModel.GetCurrentGtLikelihood(child.MaxCoverage,
                            child.GetAlleleCounts(setPosition, segmentPosition, segmentsSet), childGtState);
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

            private static void AssignMCC(PedigreeMember pedigreeMember, int setPosition, int segmentPosition, SegmentsSet segmentsSet,
                Dictionary<int, List<Genotype>> genotypes, Genotype gtStates, int copyNumber)
            {
                const int diploidCopyNumber = 2;
                const int haploidCopyNumber = 1;
                if (copyNumber > diploidCopyNumber)
                {

                    pedigreeMember.SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].MajorChromosomeCount =
                        Math.Max(gtStates.CountsA, gtStates.CountsB);
                    int? selectedGtState = genotypes[copyNumber].IndexOf(gtStates);
                    pedigreeMember.SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].MajorChromosomeCountScore =
                        pedigreeMember.CnModel.GetGtLikelihoodScore(pedigreeMember.GetAlleleCounts(setPosition, segmentPosition, segmentsSet),
                            genotypes[copyNumber], ref selectedGtState, pedigreeMember.MaxCoverage);
                }
                else
                {
                    pedigreeMember.SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].MajorChromosomeCount = copyNumber == diploidCopyNumber
                        ? haploidCopyNumber
                        : copyNumber;
                    pedigreeMember.SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].MajorChromosomeCountScore = null;
                }
            }

            private static double GetCurrentGtLikelihood(PedigreeMember pedigreeMember, int setPosition, int segmentPosition, SegmentsSet segmentsSet,
                Genotype gtStates)
            {
                double currentLikelihood = pedigreeMember.CnModel.GetCurrentGtLikelihood(pedigreeMember.MaxCoverage,
                    pedigreeMember.GetAlleleCounts(setPosition, segmentPosition, segmentsSet), gtStates);
                return currentLikelihood;
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
            public double[][] AssignCopyNumberNoPedigreeInfo(LinkedList<PedigreeMember> samples, int setPosition,
                int segmentPosition, SegmentsSet segmentsSet, List<List<int>> copyNumberCombinations)
            {
                const int defaultCn = 2;
                const double maxCoverageMultiplier = 3.0;

                double maximalLikelihood = 0;
                foreach (PedigreeMember sample in samples)
                    sample.SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].CopyNumber = defaultCn;
                int nCopies = CallerParameters.MaximumCopyNumber;
                var names = samples.Select(x => x.Name).ToList();
                var totalLikelihoods = new List<double>();
                foreach (var copyNumberCombination in copyNumberCombinations)
                {
                    double totalLikelihood = 0;
                    foreach (PedigreeMember sample in samples)
                    {
                        maximalLikelihood = 0;
                        foreach (var copyNumber in copyNumberCombination)
                        {
                            var currentLikelihood =
                                sample.CnModel.GetCnLikelihood(
                                    Math.Min(sample.GetCoverage(setPosition, segmentPosition, segmentsSet, CallerParameters.NumberOfTrimmedBins),
                                        sample.MeanCoverage * maxCoverageMultiplier))[copyNumber];
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

                var density = new double[samples.Count][];
                // no need to iterate over multiple genotypes for n (samples.Count) = 1
                if (samples.Count == 1)
                {
                    samples.First().SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].CopyNumber =
                        copyNumberCombinations[totalLikelihoods.IndexOf(totalLikelihoods.Max())].First();
                    density[0] = totalLikelihoods.ToArray();
                    return density;
                }

                var bestcopyNumberCombination = copyNumberCombinations[totalLikelihoods.IndexOf(totalLikelihoods.Max())];
                int counter = 0;
                foreach (PedigreeMember sample in samples)
                {
                    maximalLikelihood = 0;
                    density[counter] = new double[nCopies];
                    counter++;
                    foreach (var copyNumber in bestcopyNumberCombination)
                    {
                        double currentLikelihood =
                            sample.CnModel.GetCnLikelihood(
                                Math.Min(sample.GetCoverage(setPosition, segmentPosition, segmentsSet, CallerParameters.NumberOfTrimmedBins),
                                    sample.MeanCoverage * maxCoverageMultiplier))[copyNumber];
                        currentLikelihood = Double.IsNaN(currentLikelihood) || Double.IsInfinity(currentLikelihood)
                            ? 0
                            : currentLikelihood;
                        if (currentLikelihood > maximalLikelihood)
                        {
                            maximalLikelihood = currentLikelihood;
                            sample.SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].CopyNumber = copyNumber;
                            density[names.FindIndex(name => name == sample.Name)][copyNumber] = maximalLikelihood;
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
            public void AssignMccNoPedigreeInfo(LinkedList<PedigreeMember> samples, int setPosition,
                int segmentPosition, SegmentsSet segmentsSet, Dictionary<int, List<Genotype>> genotypes)
            {
                foreach (PedigreeMember sample in samples)
                {
                    int copyNumber = sample.SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].CopyNumber;
                    if (copyNumber > 2)
                    {
                        sample.SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].MajorChromosomeCount = copyNumber == 2 ? 1 : copyNumber;
                        return;
                    }
                    var genotypeset = genotypes[copyNumber];
                    int? selectedGtState = null;
                    double gqscore = sample.CnModel.GetGtLikelihoodScore(sample.GetAlleleCounts(setPosition, segmentPosition, segmentsSet),
                        genotypeset, ref selectedGtState, sample.MaxCoverage);
                    sample.SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].MajorChromosomeCountScore = gqscore;
                    if (selectedGtState.HasValue)
                        sample.SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].MajorChromosomeCount =
                            Math.Max(genotypeset[selectedGtState.Value].CountsA,
                                genotypeset[selectedGtState.Value].CountsB);
                }
            }				



            private static void InitializeCn(int setPosition, int segmentPosition,
                SegmentsSet segmentsSet, List<PedigreeMember> parents, List<PedigreeMember> children)
            {
                const int defaultCn = 2;
                parents.ForEach(sample=> sample.SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition]
                    .CopyNumber = defaultCn);
                children.ForEach(sample => sample.SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition]
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
                    return new List<int> {0};

                if (copyNumber == 1)
                    return new List<int> {0, 1};

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
                            partialGenotypes.Concat(new List<Genotype> {genotype}).ToList());
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


            public Dictionary<string, PedigreeMember.Kinship> ReadPedigreeFile(string pedigreeFile)
            {
                Dictionary<string, PedigreeMember.Kinship> kinships = new Dictionary<string, PedigreeMember.Kinship>();
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
                            kinships.Add(fields[1], PedigreeMember.Kinship.Parent);
                        else if (proband == "affected")
                            kinships.Add(fields[1], PedigreeMember.Kinship.Proband);
                        else
                            kinships.Add(fields[1], PedigreeMember.Kinship.Offspring);
                    }
                }
                return kinships;
            }

            public List<double> GetSingleSampleQualityScores(CopyNumberDistribution density, List<int> cnStates,
                List<string> sampleNames)
            {
                var singleSampleQualityScores = new List<double>();
                if (density.Count != cnStates.Count)
                    throw new ArgumentException("Size of CopyNumberDistribution should be equal to number of CN states");
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

            public double GetDeNovoQualityScore(List<PedigreeMember> parents, CopyNumberDistribution density,
                string sampleName, int sampleValue, double sampleProbability)
            {
                int nSamples = density.Count;
                const int diploidState = 2;
                var probandMarginalProbabilities = density.GetMarginalProbability(nSamples,
                    CallerParameters.MaximumCopyNumber, sampleName);
                double normalization = probandMarginalProbabilities[sampleValue] +
                                       probandMarginalProbabilities[diploidState];
                double probandMarginalAlt = probandMarginalProbabilities[sampleValue] / normalization;
                //density.SetConditionalProbability(density.Count, MaximumCopyNumber, sampleName, sampleValue, probandMarginalProbabilities[sampleValue]);

                var parentNames = parents.Select(x => x.Name).ToList();
                var firstParentMarginalProbabilities = density.GetMarginalProbability(nSamples,
                    CallerParameters.MaximumCopyNumber, parentNames.First());
                var secondParentMarginalProbabilities = density.GetMarginalProbability(nSamples,
                    CallerParameters.MaximumCopyNumber, parentNames.Last());
                normalization = firstParentMarginalProbabilities[sampleValue] +
                                firstParentMarginalProbabilities[diploidState];
                double firstParentMarginalAlt =
                    Math.Min(Math.Max(firstParentMarginalProbabilities[sampleValue] / normalization, 0.001), 0.999);
                normalization = secondParentMarginalProbabilities[sampleValue] +
                                secondParentMarginalProbabilities[diploidState];
                double secondParentMarginalAlt =
                    Math.Min(Math.Max(secondParentMarginalProbabilities[sampleValue] / normalization, 0.001), 0.999);

                normalization = (1 - firstParentMarginalAlt) * secondParentMarginalAlt +
                                firstParentMarginalAlt * secondParentMarginalAlt +
                                (1 - firstParentMarginalAlt) * (1 - firstParentMarginalAlt) +
                                (1 - secondParentMarginalAlt) * firstParentMarginalAlt;
                double diploidProbability = (1 - firstParentMarginalAlt) * (1 - secondParentMarginalAlt) / normalization;
                double denovoProbability = diploidProbability * probandMarginalAlt;
                double qscore = -10.0 * Math.Log10(1 - denovoProbability);
                return qscore;
            }

            public double GetConditionalDeNovoQualityScore(CopyNumberDistribution density, int probandIndex,
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

