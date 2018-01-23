using System.Collections.Generic;
using CanvasCommon;
using Illumina.Common.FileSystem;

namespace CanvasPedigreeCaller.Visualization
{
    public interface ICopyNumberBedGraphWriter
    {
        void Write(IReadOnlyList<CanvasSegment> segments, PloidyInfo ploidyInfo, IFileLocation location);
    }
}