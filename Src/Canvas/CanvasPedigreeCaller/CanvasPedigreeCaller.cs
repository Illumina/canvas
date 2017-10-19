using CanvasCommon;
using Combinatorics.Collections;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Isas.Framework.DataTypes;
using Isas.Framework.Logging;
using Isas.SequencingFiles;
using Isas.SequencingFiles.Vcf;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using Isas.Framework.DataTypes.Maps;

namespace CanvasPedigreeCaller
{
    class CanvasPedigreeCaller
    {
        #region Members
        public const int DefaultQualityFilterThreshold = 7;
        public const int DefaultDeNovoQualityFilterThreshold = 20;
        private readonly int _qualityFilterThreshold;
        private readonly int _deNovoQualityFilterThreshold;
        private readonly PedigreeCallerParameters _callerParameters;

        private readonly ILogger _logger;
        private readonly VariantCaller _variantCaller;
        private readonly CopyNumberLikelihoodCalculator _copyNumberLikelihoodCalculator;
        #endregion

        public CanvasPedigreeCaller(ILogger logger, int qualityFilterThreshold, int deNovoQualityFilterThreshold, PedigreeCallerParameters callerParameters, CopyNumberLikelihoodCalculator copyNumberLikelihoodCalculator, VariantCaller variantCaller)
        {
            _logger = logger;
            _qualityFilterThreshold = qualityFilterThreshold;
            _deNovoQualityFilterThreshold = deNovoQualityFilterThreshold;
            _callerParameters = callerParameters;
            _copyNumberLikelihoodCalculator = copyNumberLikelihoodCalculator;
            _variantCaller = variantCaller;
        }

        internal int CallVariants(List<string> variantFrequencyFiles, List<string> segmentFiles,
            string outVcfFile, string ploidyBedPath, string referenceFolder, List<string> sampleNames, string commonCNVsbedPath, string pedigreeFile)
        {
            // load files
            // initialize data structures and classes
            var fileCounter = 0;
            var samplesInfo = new SampleMap<SampleMetrics>();
            var sampleSegments = new SampleMap<Segments>();
            var copyNumberModels = new SampleMap<CopyNumberModel>();
            var variantFrequencyFilesSampleList = new SampleMap<string>();
            ISampleMap<Kinship> kinships;
            if (pedigreeFile != null)
            {
                kinships = ReadPedigreeFile(pedigreeFile);
                // In Kinship enum Proband gets the highest int value 
                var newSampleNames = kinships.OrderByDescending(x => x.Value).Select(x => x.Key.ToString()).ToList();
                var remapIndices = newSampleNames.Select(newname => sampleNames.FindIndex(name => name == newname)).ToList();
                segmentFiles = remapIndices.Select(index => segmentFiles[index]).ToList();
                variantFrequencyFiles = remapIndices.Select(index => variantFrequencyFiles[index]).ToList();
                sampleNames = newSampleNames;
            }
            else
            {
                kinships = sampleNames.Select(name => (new SampleId(name), Kinship.Other)).ToSampleMap();
            }

            foreach (string sampleName in sampleNames)
            {
                var sampleId = new SampleId(sampleName);
                var segment = Segments.ReadSegments(_logger, new FileLocation(segmentFiles[fileCounter]));
                segment.AddAlleles(CanvasIO.ReadFrequenciesWrapper(_logger, new FileLocation(variantFrequencyFiles[fileCounter]), segment.IntervalsByChromosome));
                sampleSegments.Add(sampleId, segment);
                var sampleInfo = SampleMetrics.GetSampleInfo(segment, ploidyBedPath, _callerParameters.NumberOfTrimmedBins, sampleId);
                var copyNumberModel = new CopyNumberModel(_callerParameters.MaximumCopyNumber, sampleInfo.MeanMafCoverage, sampleInfo.MeanCoverage, sampleInfo.MaxCoverage);
                samplesInfo.Add(sampleId, sampleInfo);
                copyNumberModels.Add(sampleId, copyNumberModel);
                variantFrequencyFilesSampleList.Add(sampleId, variantFrequencyFiles[fileCounter]);
                fileCounter++;
            }
            var segmentSetsFromCommonCnvs = CreateSegmentSetsFromCommonCnvs(variantFrequencyFilesSampleList,
                _callerParameters.DefaultReadCountsThreshold, commonCNVsbedPath, sampleSegments);

            var segmentsForVariantCalling = GetHighestLikelihoodSegments(segmentSetsFromCommonCnvs, samplesInfo, copyNumberModels).ToList();
            PedigreeInfo pedigreeInfo = null;
            if (kinships.Values.Any(kin => kin == Kinship.Proband))
                pedigreeInfo = PedigreeInfo.GetPedigreeInfo(kinships, _callerParameters);
            Parallel.ForEach(
                segmentsForVariantCalling,
                new ParallelOptions
                {
                    MaxDegreeOfParallelism = Math.Min(Environment.ProcessorCount, _callerParameters.MaxCoreNumber)
                },
                segments => _variantCaller.CallVariant(segments, samplesInfo, copyNumberModels, pedigreeInfo)
            );
            var variantCalledSegments = new SampleMap<List<CanvasSegment>>();
            foreach (var key in samplesInfo.SampleIds)
                variantCalledSegments.Add(key, segmentsForVariantCalling.Select(segment => segment[key]).ToList());

            var mergedVariantCalledSegments = MergeSegments(variantCalledSegments, _callerParameters.MinimumCallSize);
            var outputFolder = new FileLocation(outVcfFile).Directory;
            foreach (var sampleId in samplesInfo.SampleIds)
            {
                var coverageOutputPath = SingleSampleCallset.GetCoverageAndVariantFrequencyOutput(outputFolder,
                    sampleId.ToString());
                CanvasSegment.WriteCoveragePlotData(mergedVariantCalledSegments[sampleId], samplesInfo[sampleId].MeanCoverage,
                    samplesInfo[sampleId].Ploidy, coverageOutputPath, referenceFolder);
            }
            bool isPedigreeInfoSupplied = pedigreeInfo != null;
            var denovoQualityThreshold = isPedigreeInfoSupplied ? (int?)_deNovoQualityFilterThreshold : null;
            var ploidies = samplesInfo.Select(info => info.Value.Ploidy).ToList();
            var diploidCoverage = samplesInfo.Select(info => info.Value.MeanCoverage).ToList();
            var names = samplesInfo.SampleIds.Select(id => id.ToString()).ToList();
            CanvasSegmentWriter.WriteMultiSampleSegments(outVcfFile, mergedVariantCalledSegments, diploidCoverage, referenceFolder, names,
                null, ploidies, _qualityFilterThreshold, isPedigreeInfoSupplied, denovoQualityThreshold);

            outputFolder = new FileLocation(outVcfFile).Directory;
            foreach (var sampleId in samplesInfo.SampleIds)
            {
                var outputVcfPath = SingleSampleCallset.GetSingleSamplePedigreeVcfOutput(outputFolder, sampleId.ToString());
                CanvasSegmentWriter.WriteSegments(outputVcfPath.FullName, mergedVariantCalledSegments[sampleId],
                    samplesInfo[sampleId].MeanCoverage, referenceFolder, sampleId.ToString(), null,
                    samplesInfo[sampleId].Ploidy, _qualityFilterThreshold, isPedigreeInfoSupplied, denovoQualityThreshold);
            }
            return 0;
        }

        private IEnumerable<ISampleMap<CanvasSegment>> GetHighestLikelihoodSegments(IEnumerable<ISampleMap<OverlappingSegmentsRegion>> segmentSetsFromCommonCnvs,
            ISampleMap<SampleMetrics> pedigreeMembersInfo, ISampleMap<CopyNumberModel> copyNumberModel)
        {
            var updatedSegmentSets = segmentSetsFromCommonCnvs
                .AsParallel()
                .AsOrdered()
                .WithDegreeOfParallelism(Math.Min(Environment.ProcessorCount, _callerParameters.MaxCoreNumber))
                .Select(segmentSet =>
                {
                    GetHighestLikelihoodSegmentsSet(segmentSet, pedigreeMembersInfo, copyNumberModel);
                    return segmentSet;
                });

            return updatedSegmentSets
                .SelectMany(sampleMap => sampleMap.SelectValues(x => x.GetSet().AsEnumerable()).Zip())
                .ToList();
        }


        private static ISampleMap<List<CanvasSegment>> MergeSegments(ISampleMap<List<CanvasSegment>> segments, int minimumCallSize)
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

            var mergedSegments = new SampleMap<List<CanvasSegment>>();
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
        private IEnumerable<ISampleMap<OverlappingSegmentsRegion>> CreateSegmentSetsFromCommonCnvs(ISampleMap<string> variantFrequencyFiles,
            int defaultAlleleCountThreshold, string commonCNVsbedPath, ISampleMap<Segments> sampleSegments)
        {
            if (commonCNVsbedPath == null)
            {
                var defaultSampleRegions = sampleSegments
                    .SelectValues(segments => segments.AllSegments.Select(segment => new OverlappingSegmentsRegion(segment)).ToList());
                return GetOverlappingSegmentsRegionSampleLists(defaultSampleRegions);
            }

            var commonRegions = ReadCommonRegions(commonCNVsbedPath);
            var chromosomes = sampleSegments.Values.First().GetChromosomes();
            if (IsIdenticalChromosomeNames(commonRegions, chromosomes))
                throw new ArgumentException(
                    $"Chromosome names in a common CNVs bed file {commonCNVsbedPath} does not match the genome reference");

            var segmentIntervalsByChromosome = new Dictionary<string, List<BedInterval>>();
            var genomicBinsByChromosome = new Dictionary<string, IReadOnlyList<SampleGenomicBin>>();

            Parallel.ForEach(
                chromosomes,
                chr =>
                {
                    genomicBinsByChromosome[chr] = sampleSegments.Values.First().GetGenomicBinsForChromosome(chr);
                    segmentIntervalsByChromosome[chr] =
                        CanvasSegment.RemapGenomicToBinCoordinates(commonRegions[chr], genomicBinsByChromosome[chr]);
                });

            var sampleRegions = new SampleMap<List<OverlappingSegmentsRegion>>();
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

        private static IEnumerable<ISampleMap<OverlappingSegmentsRegion>> GetOverlappingSegmentsRegionSampleLists(ISampleMap<List<OverlappingSegmentsRegion>> sampleRegions)
        {
            return sampleRegions.Zip();
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

        public static int AggregateVariantCoverage(ref List<CanvasSegment> segments)
        {
            var variantCoverage = segments.SelectMany(segment => segment.Balleles.TotalCoverage).ToList();
            return variantCoverage.Any() ? Utilities.Median(variantCoverage) : 0;
        }

        /// <summary>
        /// Identify variant with the highest likelihood at a given setPosition and assign relevant scores
        /// </summary>
        private void GetHighestLikelihoodSegmentsSet(ISampleMap<OverlappingSegmentsRegion> canvasSegmentsSet, ISampleMap<SampleMetrics> pedigreeMembersInfo,
            ISampleMap<CopyNumberModel> model)
        {
            SegmentsSet segmentSet;

            if (canvasSegmentsSet.Values.First().SetA == null)
                segmentSet = SegmentsSet.SetB;
            else if (canvasSegmentsSet.Values.First().SetB == null)
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

        private double GetSegmentSetLikelihood(ISampleMap<OverlappingSegmentsRegion> canvasSegmentsSet, ISampleMap<SampleMetrics> samplesInfo,
            ISampleMap<CopyNumberModel> copyNumberModel, SegmentsSet segmentsSet)
        {
            double segmentSetLikelihood = 0;
            foreach (var sampleId in canvasSegmentsSet.SampleIds)
                canvasSegmentsSet[sampleId].SetSet(segmentsSet);

            var canvasSegments = new List<ISampleMap<CanvasSegment>>();
            int nSegments = canvasSegmentsSet.First().Value.GetSet().Count;
            for (var canvasSegmentIndex = 0; canvasSegmentIndex < nSegments; canvasSegmentIndex++)
            {
                var canvasSegment = new SampleMap<CanvasSegment>();
                foreach (var id in canvasSegmentsSet.SampleIds)
                    canvasSegment.Add(id, canvasSegmentsSet[id].GetSet()[canvasSegmentIndex]);
                canvasSegments.Add(canvasSegment);
            }
            foreach (var canvasSegment in canvasSegments)
            {
                var copyNumbersLikelihoods = _copyNumberLikelihoodCalculator.GetCopyNumbersLikelihoods(canvasSegment, samplesInfo, copyNumberModel);
                GetCopyNumbersNoPedigreeInfo(canvasSegment, copyNumbersLikelihoods);
                segmentSetLikelihood += copyNumbersLikelihoods.MaximalLikelihood;
            }

            segmentSetLikelihood /= nSegments;
            return segmentSetLikelihood;
        }

        /// <summary>
        /// Calculates maximal likelihood for copy numbers. Updated CanvasSegment CopyNumber only. 
        /// </summary>
        public static ISampleMap<int> GetCopyNumbersNoPedigreeInfo(ISampleMap<CanvasSegment> segments, CopyNumbersLikelihoods copyNumbersLikelihoods)
        {
            var sampleCopyNumbers = new SampleMap<int>();
            double maximalLikelihood = 1;
            foreach (var sampleId in segments.SampleIds)
            {
                var (copyNumber, maxSampleLikelihood) = MaxBy(copyNumbersLikelihoods.SingleSampleLikelihoods[sampleId], x => x.Value);
                maximalLikelihood *= maxSampleLikelihood;
                sampleCopyNumbers.Add(sampleId, copyNumber);
            }
            copyNumbersLikelihoods.MaximalLikelihood = maximalLikelihood;
            return sampleCopyNumbers;
        }

        public static T MaxBy<T>(IEnumerable<T> items, Func<T, double> transform)
        {
            bool any = false;
            var maxValue = double.MinValue;
            var itemWithMaxValue = default(T);
            foreach (var item in items)
            {
                any = true;
                var value = transform(item);
                if (value <= maxValue) continue;
                itemWithMaxValue = item;
                maxValue = value;
            }
            if (!any)
                throw new IlluminaException("Cannot get max value from an empty enumerable");
            return itemWithMaxValue;
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

        public double GetTransitionProbability(int gt1Parent, int gt2Parent, int gt1Offspring, int gt2Offspring)
        {
            if (gt1Parent == gt1Offspring || gt1Parent == gt2Offspring ||
                gt2Parent == gt1Offspring || gt2Parent == gt2Offspring)
                return 0.5;
            return _callerParameters.DeNovoRate;
        }

        public enum Kinship
        {
            Other = 0,
            Parent = 1,
            Proband = 2
        }

        public static ISampleMap<Kinship> ReadPedigreeFile(string pedigreeFile)
        {
            var kinships = new SampleMap<Kinship>();
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
    }
}