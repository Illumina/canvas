using System.Collections.Generic;
using CanvasCommon;
using Isas.SequencingFiles.Bed;

namespace CanvasPedigreeCaller.Visualization
{
    public interface IBedGraphCalculator
    {
        IEnumerable<BedGraphEntry> Calculate(IReadOnlyList<CanvasSegment> segments);
    }
}