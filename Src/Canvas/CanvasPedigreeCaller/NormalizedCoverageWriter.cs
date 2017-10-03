using System;
using System.Collections.Generic;
using System.Linq;
using CanvasCommon;
using Illumina.Common.FileSystem;
using Illumina.Common.MathUtilities;
using Isas.SequencingFiles.Bed;

namespace CanvasPedigreeCaller
{
    public class NormalizedCoverageWriter
    {
        public static void Write(IReadOnlyList<CanvasSegment> segments, IFileLocation coverageFile)
        {
            var normalizationFactor = ComputeNormalizationFactor(GetSegmentsForNormalizationEstimation(segments));
            var begGraphEntries = segments.SelectMany<CanvasSegment, BedGraphEntry>(segment => GetBedGraphEntries(segment, normalizationFactor));
            BedGraphWriter.WriteLines(coverageFile, begGraphEntries);
        }

        private static IEnumerable<BedGraphEntry> GetBedGraphEntries(CanvasSegment segment, double normalizationFactor)
        {
            return segment.GenomicBins.Select<SampleGenomicBin, BedGraphEntry>(bin => GetBedGraphEntry(bin, normalizationFactor));
        }

        private static BedGraphEntry GetBedGraphEntry(SampleGenomicBin bin, double normalizationFactor)
        {
            var normalizedCoverage = bin.Count * normalizationFactor;
            return new BedGraphEntry(bin.GenomicBin.Chromosome, bin.GenomicBin.Interval, (decimal)normalizedCoverage);
        }

        private static IEnumerable<CanvasSegment> GetSegmentsForNormalizationEstimation(IReadOnlyList<CanvasSegment> segments)
        {
            bool PassFilter(CanvasSegment segment) => segment.Filter == "PASS";
            if (segments.Any(PassFilter))
                return segments.Where(PassFilter);
            return segments;
        }

        private static double ComputeNormalizationFactor(IEnumerable<CanvasSegment> segments)
        {
            var weightedNormalizationFactors = ParallelEnumerable.Select<CanvasSegment, ValueTuple<double, double>>(segments
                    .AsParallel()
                    .AsOrdered()
                    .Where(segment => segment.CopyNumber != 0), segment => (GetNormalizationFactor(segment), GetWeight(segment)));
            return WeightedMedian.Median(weightedNormalizationFactors);
        }

        private static double GetWeight(CanvasSegment segment)
        {
            return segment.BinCount;
        }

        private static double GetNormalizationFactor(CanvasSegment segment)
        {
            return segment.CopyNumber / WeightedMedian.Median(segment.Counts.Select(count => ((double)count, 1d)));
        }
    }
}