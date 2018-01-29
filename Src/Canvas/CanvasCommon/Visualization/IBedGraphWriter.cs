using System.Collections.Generic;
using Illumina.Common.FileSystem;
using Isas.SequencingFiles.Bed;

namespace CanvasCommon.Visualization
{
    public interface IBedGraphWriter
    {
        void Write(IEnumerable<BedGraphEntry> bedGraphEntries, IFileLocation location);
    }
}