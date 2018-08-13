using System.Collections.Generic;
using System.Linq;
using CanvasCommon;
using Isas.SequencingFiles;
using Isas.SequencingFiles.Bed;
using MathNet.Numerics.Statistics;

namespace CanvasPedigreeCaller.Visualization
{
    public class NormalizedSegmentsCoverageCalculator : BaseNormalizedCoverageCalculator
    {
        protected override IEnumerable<BedGraphEntry> GetBedGraphEntries(CanvasSegment segment, double normalizationFactor)
        {
            return new List<BedGraphEntry>(){GetBedGraphEntry(segment, normalizationFactor)};
        }

        private static BedGraphEntry GetBedGraphEntry(CanvasSegment segment, double normalizationFactor)
        {
            var segmentBins = segment.GenomicBins;
            var medianBinCount = segmentBins.Select(b => b.Count).Median();

            var normalizedCoverage = (decimal)(medianBinCount * normalizationFactor);
            return new BedGraphEntry(segment.Chr, new BedInterval(segmentBins.Min(x=>x.Start), segmentBins.Max(x=>x.Stop)), normalizedCoverage);
        }

    }
}