using System.Collections.Generic;
using CanvasCommon;
using CanvasCommon.Visualization;
using Illumina.Common.FileSystem;
using Isas.Framework.DataTypes;
using Isas.SequencingFiles;
using Isas.SequencingFiles.Bed;

namespace CanvasPedigreeCaller.Visualization
{
    public interface ICoverageBigWigWriter
    {
        IFileLocation Write(IReadOnlyList<CanvasSegment> segments, IDirectoryLocation output,
            double normalizationFactor);
    }
}