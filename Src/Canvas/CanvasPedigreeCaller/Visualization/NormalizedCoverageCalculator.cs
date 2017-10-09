using System.Collections.Generic;
using System.Linq;
using CanvasCommon;
using Illumina.Common.MathUtilities;
using Isas.SequencingFiles.Bed;

namespace CanvasPedigreeCaller.Visualization
{
    public class NormalizedCoverageCalculator
    {
        public IEnumerable<BedGraphEntry> Calculate(IReadOnlyList<CanvasSegment> segments)
        {
            var normalizationFactor = ComputeNormalizationFactor(segments);
            return segments.SelectMany(segment => GetBedGraphEntries(segment, normalizationFactor));
        }

        private static IEnumerable<BedGraphEntry> GetBedGraphEntries(CanvasSegment segment, double normalizationFactor)
        {
            return segment.GenomicBins.Select(bin => GetBedGraphEntry(bin, normalizationFactor));
        }

        private static BedGraphEntry GetBedGraphEntry(SampleGenomicBin bin, double normalizationFactor)
        {
            var normalizedCoverage = (decimal)(bin.Count * normalizationFactor);
            return new BedGraphEntry(bin.GenomicBin.Chromosome, bin.GenomicBin.Interval, normalizedCoverage);
        }

        private static double ComputeNormalizationFactor(IReadOnlyList<CanvasSegment> segments)
        {
            var segmentsForEstimation = GetSegmentsForNormalizationEstimation(segments);
            var weightedNormalizationFactors = segmentsForEstimation
                .AsParallel()
                .Where(segment => segment.CopyNumber != 0)
                .Select(segment => (GetNormalizationFactor(segment), GetWeight(segment)));
            return WeightedMedian.Median(weightedNormalizationFactors);
        }

        private static IEnumerable<CanvasSegment> GetSegmentsForNormalizationEstimation(IReadOnlyList<CanvasSegment> segments)
        {
            bool PassFilter(CanvasSegment segment) => segment.Filter == "PASS";
            if (segments.Any(PassFilter))
                return segments.Where(PassFilter);
            return segments;
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