using System.Collections.Generic;
using CanvasCommon;
using Illumina.Common.FileSystem;

namespace CanvasPedigreeCaller.Visualization
{
    public interface ICoverageBigWigWriter
    {
        IFileLocation Write(IReadOnlyList<CanvasSegment> segments, IDirectoryLocation output);
    }
}