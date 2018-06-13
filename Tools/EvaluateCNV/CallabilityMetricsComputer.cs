using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Isas.Framework.Logging;
using Isas.Framework.Utilities;
using Isas.SequencingFiles;
using Isas.SequencingFiles.Vcf;
using Isas.StatisticsProcessing.Processors.StatsProcessor;

namespace EvaluateCNV
{
    public class CallabilityMetricsComputer
    {
        private readonly ILogger _logger;
        private readonly CallabilityCalculator _callability;
        private readonly StatsProcessor<ReferenceInterval> _bedProcessor;

        private CallabilityMetricsComputer(ILogger logger, CallabilityCalculator callability)
        {
            _logger = logger;
            _callability = callability;
            var metricsCollection = new StatCalculatorCollection<ReferenceInterval>
            {
                _callability,
            };
            _bedProcessor = new StatsProcessor<ReferenceInterval>(metricsCollection);
        }

        public CallabilityMetric CalculateMetric(Dictionary<string, List<CnvCall>> calls)
        {
            var callIntervals = calls.SelectMany(chrCalls => chrCalls.Value).Select(call =>
                new ReferenceInterval(call.Chr, new Interval(call.Start, call.End)));
            _logger.Info("Running callability calculation");
            var benchmark = new Benchmark();
            var results = _bedProcessor.Process(callIntervals);
            var callability = results.GetResult(_callability);
            _logger.Info($"Callability calculation elapsed time: {benchmark.GetElapsedTime()}");
            return callability;
        }

        public static CallabilityMetricsComputer Create(ILogger logger, GenomeMetadata genomeMetadata, IFileLocation filterBed, bool isFemale)
        {
            var filterIntervals = LoadBedRegions(filterBed, genomeMetadata);
            var nonFilterIntervals = GetIncludedIntervals(filterIntervals, genomeMetadata);
            string chrY = genomeMetadata
                .Contigs()
                .Select(contig => contig.Name)
                .SingleOrDefault(name => name == "chrY" || name == "Y");
            if (isFemale && chrY != null)
                nonFilterIntervals.Remove(chrY);
            var callabilityCalculator = new CallabilityCalculator(nonFilterIntervals);

            return new CallabilityMetricsComputer(logger, callabilityCalculator);
        }

        private static Dictionary<string, IEnumerable<Interval>> GetIncludedIntervals(Dictionary<string, List<Interval>> filterIntervals, GenomeMetadata genomeMetadata)
        {
            var contigNames = genomeMetadata
                .Contigs()
                .Where(contig => contig.Type == GenomeMetadata.SequenceType.Autosome || contig.Type == GenomeMetadata.SequenceType.Allosome)
                .Select(contig => contig.Name);
            return contigNames
                .Select(contig =>
                {
                    var filteredIntervals = filterIntervals.ContainsKey(contig)
                        ? filterIntervals[contig]
                        : new List<Interval>();
                    return (contig, GetIncludedIntervals(filteredIntervals, genomeMetadata.GetSequence(contig).Length));
                })
                .ToDictionary();
        }

        private static IEnumerable<Interval> GetIncludedIntervals(List<Interval> intervals, long totalLength)
        {
            int currentPosition = 1;
            foreach (var interval in intervals)
            {
                if (currentPosition < interval.OneBasedStart)
                {
                    yield return new Interval(currentPosition, interval.OneBasedStart - 1);
                }
                currentPosition = interval.OneBasedEnd + 1;
            }
            if (currentPosition <= totalLength)
                yield return new Interval(currentPosition, (int)totalLength);
        }

        private static Dictionary<string, List<Interval>> LoadBedRegions(IFileLocation bedFile,
            GenomeMetadata genomeMetadata)
        {
            return File.ReadAllLines(bedFile.FullName)
                .Select(line =>
                {
                    var bits = line.Split('\t');
                    return (Chromosome: bits[0], Interval: new Interval(int.Parse(bits[1]) + 1, int.Parse(bits[2])));
                })
                .GroupByAdjacent(bedEntry => bedEntry.Chromosome)
                .Where(kvp => genomeMetadata.GetSequence(kvp.Key) != null)
                .ToDictionary(
                    chromosomeElements => chromosomeElements.Key,
                    chromosomeElements => GetValidatedIntervals(chromosomeElements.Value,
                        genomeMetadata.GetSequence(chromosomeElements.Key).Length).ToList());
        }

        private static IEnumerable<Interval> GetValidatedIntervals(
            IEnumerable<(string Chromosome, Interval Interval)> intervals, long chrLength)
        {
            Interval currentInterval = null;
            foreach (var (chr, interval) in intervals)
            {
                if (interval.OneBasedStart < 1 || interval.OneBasedEnd > chrLength)
                    throw new ArgumentException($"Invalid interval for {chr} with length {chrLength}: {interval}");
                if (currentInterval != null && interval.OneBasedStart < currentInterval.OneBasedStart)
                    throw new ArgumentException(
                        $"Unsorted intervals for chr {chr}: '{currentInterval}' followed by '{interval}'");
                yield return interval;
                currentInterval = interval;
            }
        }
    }
}