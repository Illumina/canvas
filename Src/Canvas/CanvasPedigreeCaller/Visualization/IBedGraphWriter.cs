using System.Collections.Generic;
using CanvasCommon;
using Illumina.Common.FileSystem;
using Isas.SequencingFiles.Bed;

namespace CanvasPedigreeCaller.Visualization
{
    public interface IBedGraphWriter
    {
        void Write(IEnumerable<BedGraphEntry> bedGraphEntries, IFileLocation location);
    }
}