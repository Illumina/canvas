using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Isas.Framework.Logging;
using Isas.SequencingFiles;
using Isas.SequencingFiles.Vcf;


namespace CanvasCommon
{
    public class Segments
    {
        private readonly OrderedDictionary<string, List<CanvasSegment>> _segments;

        public readonly IReadOnlyList<CanvasSegment> AllSegments;

        private Segments(OrderedDictionary<string, List<CanvasSegment>> segments, List<CanvasSegment> allSegments)
        {
            _segments = segments;
            AllSegments = new ReadOnlyCollection<CanvasSegment>(allSegments);
        }
        public ICollection<string> GetChromosomes()
        {
            return _segments.Keys;
        }
        public IReadOnlyList<CanvasSegment> GetSegmentsForChromosome(string chromosome)
        {
            return _segments[chromosome];
        }

        public static Segments ReadSegments(ILogger logger, IFileLocation partitionedFile)
        {
            using (var reader = new GzipOrTextReader(partitionedFile.FullName))
            {
                logger.Info($"Read segments from {partitionedFile}");
                var segments = ReadSegments(reader);
                logger.Info($"Loaded {segments.AllSegments.Count} segments");
                return segments;
            }
        }

        public static Segments ReadSegments(GzipOrTextReader partitionedFileReader)
        {
            var rows = partitionedFileReader.ReadLines().Select(line => line.Split('\t'));
            var rowsByChr = rows.GroupByAdjacent(GetChromosome);
            var segments = rowsByChr.ToOrderedDictionary(kvp => kvp.Key, kvp => CreateSegments(kvp.Value));

            var allSegments = segments.SelectMany(kvp => kvp.Value).ToList();
            return new Segments(segments, allSegments);
        }

        private static List<CanvasSegment> CreateSegments(IEnumerable<string[]> rows)
        {
            var groupedBins = rows.GroupByAdjacent(GetSegmentId)
                .Select(kvp => kvp.Value.Select(CreateBin).ToList())
                .ToList();

            return groupedBins.Select((bins, index) =>
            {
                var startBin = bins.First().GenomicBin;
                var previousBin = index == 0 ? null : groupedBins[index - 1].Last().GenomicBin;
                var startConfidenceInterval = GetStartConfidenceInterval(startBin, previousBin);

                var endBin = bins.Last().GenomicBin;
                var nextBin = index == groupedBins.Count - 1 ? null : groupedBins[index + 1].First().GenomicBin;
                var endConfidenceInterval = GetEndConfidenceInterval(endBin, nextBin);

                return CreateSegment(bins, startConfidenceInterval, endConfidenceInterval);
            }).ToList();
        }

        private static Tuple<int, int> GetStartConfidenceInterval(GenomicBin thisBin, GenomicBin previousBin)
        {
            if (previousBin == null || previousBin.Interval.End != thisBin.Interval.Start)
                return GetConfidenceInterval(thisBin);
            return Tuple.Create(-GetHalfLength(previousBin), GetHalfLength(thisBin));
        }

        private static Tuple<int, int> GetEndConfidenceInterval(GenomicBin thisBin, GenomicBin nextBin)
        {
            if (nextBin == null || thisBin.Interval.End != nextBin.Interval.Start)
                return GetConfidenceInterval(thisBin);
            return Tuple.Create(-GetHalfLength(thisBin), GetHalfLength(nextBin));
        }

        private static int GetHalfLength(GenomicBin bin)
        {
            return (int)Math.Round(bin.Interval.Length / 2.0, MidpointRounding.AwayFromZero);
        }

        private static Tuple<int, int> GetConfidenceInterval(GenomicBin bin)
        {
            var halfLength = GetHalfLength(bin);
            return Tuple.Create(-halfLength, halfLength);
        }

        private static CanvasSegment CreateSegment(List<SampleGenomicBin> bins, Tuple<int, int> startConfidenceInterval, Tuple<int, int> endConfidenceInterval)
        {
            var segment = new CanvasSegment(bins.First().GenomicBin.Chromosome, bins.First().GenomicBin.Interval.Start, bins.Last().GenomicBin.Interval.End,
                bins)
            {
                StartConfidenceInterval = startConfidenceInterval,
                EndConfidenceInterval = endConfidenceInterval
            };
            return segment;
        }

        private static SampleGenomicBin CreateBin(string[] row)
        {
            var chr = GetChromosome(row);
            int start = GetStart(row);
            int end = GetEnd(row);
            var count = float.Parse(row[3]);
            return new SampleGenomicBin(chr, start, end, count);
        }

        private static int GetEnd(string[] row)
        {
            return int.Parse(row[2]);
        }

        private static int GetStart(string[] row)
        {
            return int.Parse(row[1]);
        }

        private static string GetChromosome(string[] row)
        {
            return row[0];
        }

        private static string GetSegmentId(string[] row)
        {
            return row[4];
        }

        public void AddAlleles(Dictionary<string, List<Balleles>> allelesByChromosome)
        {
            foreach (string chr in GetChromosomes())
            {
                var counter = 0;
                foreach (var segment in GetSegmentsForChromosome(chr))
                {
                    segment.Balleles.Add(allelesByChromosome[chr][counter]);
                    counter++;
                }
            }
        }
        public  Dictionary<string, List<BedInterval>> GetIntervalsByChromosome()
        {
            var intervalsByChromosome = new Dictionary<string, List<BedInterval>>();
            foreach (string chr in GetChromosomes())
            {
                intervalsByChromosome[chr] = new List<BedInterval>();
                foreach (var canvasSegment in GetSegmentsForChromosome(chr))
                {
                    intervalsByChromosome[chr].Add(new BedInterval(canvasSegment.Begin, canvasSegment.End));
                }
            }
            return intervalsByChromosome;
        }
    }
}
