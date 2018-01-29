using System.Collections.Generic;
using CanvasCommon;
using Isas.Framework.DataTypes;

namespace CanvasPedigreeCaller.Visualization
{
    public interface ICopyNumberBedGraphWriter
    {
        void Write(IReadOnlyList<CanvasSegment> segments, PloidyInfo ploidyInfo, BgzfFile location);
    }
}