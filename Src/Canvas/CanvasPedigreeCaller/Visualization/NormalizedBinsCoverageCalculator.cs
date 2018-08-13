using System.Collections.Generic;
using System.Linq;
using CanvasCommon;
using Isas.SequencingFiles.Bed;

namespace CanvasPedigreeCaller.Visualization
{
    public class NormalizedBinsCoverageCalculator : BaseNormalizedCoverageCalculator
    {
        protected override IEnumerable<BedGraphEntry> GetBedGraphEntries(CanvasSegment segment, double normalizationFactor)
        {
            return segment.GenomicBins.Select(bin => GetBedGraphEntry(bin, normalizationFactor));
        }

        private static BedGraphEntry GetBedGraphEntry(SampleGenomicBin bin, double normalizationFactor)
        {
            var normalizedCoverage = (decimal)(bin.Count * normalizationFactor);
            return new BedGraphEntry(bin.GenomicBin.Chromosome, bin.GenomicBin.Interval, normalizedCoverage);
        }

    }
}