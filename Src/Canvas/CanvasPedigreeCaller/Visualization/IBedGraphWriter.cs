using System.Collections.Generic;
using CanvasCommon;
using Illumina.Common.FileSystem;

namespace CanvasPedigreeCaller.Visualization
{
    public interface IBedGraphWriter
    {
        void Write(IReadOnlyList<CanvasSegment> segments, IFileLocation location);
    }
}