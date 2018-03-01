using System.Collections.Generic;
using Illumina.Common.FileSystem;
using Isas.SequencingFiles;
using Isas.SequencingFiles.Bed;

namespace CanvasCommon.Visualization
{
    public class BedGraphWriterFacade : IBedGraphWriter
    {
        public void Write(IEnumerable<BedGraphEntry> bedGraphEntries, IFileLocation location)
        {
            using (var bgzipWriter = new BgzipOrStreamWriter(location.FullName))
            {
                var writer = new BedGraphWriter(bgzipWriter);
                writer.WriteLines(bedGraphEntries);
            }
        }
    }
}