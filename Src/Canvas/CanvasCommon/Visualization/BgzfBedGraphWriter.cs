using System.Collections.Generic;
using Isas.ClassicBioinfoTools.Tabix;
using Isas.Framework.DataTypes;
using Isas.SequencingFiles.Bed;

namespace CanvasCommon.Visualization
{
    public class BgzfBedGraphWriter : IBgzfBedGraphWriter
    {
        private readonly IBedGraphWriter _bedGraphWriter;
        private readonly ITabixWrapper _tabixWrapper;

        public BgzfBedGraphWriter(IBedGraphWriter bedGraphWriter, ITabixWrapper tabixWrapper)
        {
            _bedGraphWriter = bedGraphWriter;
            _tabixWrapper = tabixWrapper;
        }
        public void Write(IEnumerable<BedGraphEntry> bedGraphEntries, BgzfFile location)
        {
            _bedGraphWriter.Write(bedGraphEntries, location.FileLocation);
            _tabixWrapper.BuildTabixIndex(location, TabixFileType.Bed);
        }
    }
}