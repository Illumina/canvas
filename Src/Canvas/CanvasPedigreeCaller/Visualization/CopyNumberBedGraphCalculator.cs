using System.Collections.Generic;
using System.Linq;
using CanvasCommon;
using Isas.SequencingFiles;
using Isas.SequencingFiles.Bed;

namespace CanvasPedigreeCaller.Visualization
{
    public class CopyNumberBedGraphCalculator
    {
        public IEnumerable<BedGraphEntry> Calculate(IReadOnlyList<CanvasSegment> segments, PloidyInfo ploidyInfo)
        {
            var variantSegments = segments.Where(segment => IsPassVariant(segment, ploidyInfo));
            return variantSegments.Select(GetCopyNumberEntry);
        }

        private bool IsPassVariant(CanvasSegment segment, PloidyInfo ploidyInfo)
        {
            if (!segment.Filter.IsPass) return false;
            var referenceCopyNumber = ploidyInfo?.GetReferenceCopyNumber(segment) ?? 2;
            if (segment.CopyNumber != referenceCopyNumber) return true;
            if (segment.CopyNumber == segment.MajorChromosomeCount) return true; //LOH
            return false;
        }

        private BedGraphEntry GetCopyNumberEntry(CanvasSegment segment)
        {
            var interval = new BedInterval(segment.Begin, segment.End);
            return new BedGraphEntry(segment.Chr, interval, segment.CopyNumber);
        }
    }
}