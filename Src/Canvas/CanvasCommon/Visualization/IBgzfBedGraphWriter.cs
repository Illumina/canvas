using System.Collections.Generic;
using Isas.Framework.DataTypes;
using Isas.SequencingFiles.Bed;

namespace CanvasCommon.Visualization
{
    public interface IBgzfBedGraphWriter
    {
        void Write(IEnumerable<BedGraphEntry> bedGraphEntries, BgzfFile location, string header = null);
    }
}