using System;
using System.Collections.Generic;
using System.Linq;
using CanvasCommon;
using Illumina.Common;
using Illumina.Common.MathUtilities;

namespace CanvasPedigreeCaller.Visualization
{
    public static class NormalizationCalculator
    {
        public static double ComputeNormalizationFactor(IReadOnlyList<CanvasSegment> segments)
        {
            var segmentsForEstimation = Enumerable.Where<CanvasSegment>(GetSegmentsForNormalizationEstimation(segments), segment => segment.CopyNumber != 0).ToReadOnlyList();
            if (!segmentsForEstimation.Any())
            {
                // if all segments have CN=0 then use normalization factor of 0
                return 0;
            }
            var weightedNormalizationFactors = ParallelEnumerable.Select<CanvasSegment, (double, double)>(segmentsForEstimation
                .AsParallel(), segment => (GetNormalizationFactor(segment), GetWeight(segment)));

            return WeightedMedian.Median(weightedNormalizationFactors);
        }

        private static IEnumerable<CanvasSegment> GetSegmentsForNormalizationEstimation(IReadOnlyList<CanvasSegment> segments)
        {
            bool PassFilter(CanvasSegment segment) => segment.Filter.IsPass;
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