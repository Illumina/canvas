using System.Collections.Generic;
using Isas.Framework.DataTypes;
using Isas.SequencingFiles.Bed;

namespace CanvasPedigreeCaller.Visualization
{
    public interface IBgzfBedGraphWriter
    {
        void Write(IEnumerable<BedGraphEntry> bedGraphEntries, BgzfFile location);
    }
}