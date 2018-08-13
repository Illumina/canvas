using System.Collections.Generic;
using Illumina.Common.FileSystem;
using Isas.SequencingFiles;
using Isas.SequencingFiles.Bed;

namespace CanvasCommon.Visualization
{
    public class BedGraphWriterFacade : IBedGraphWriter
    {
        public void Write(IEnumerable<BedGraphEntry> bedGraphEntries, IFileLocation location, string header = null)
        {
            using (var bgzipWriter = new BgzipOrStreamWriter(location.FullName))
            {
                if (header != null)
                {
                    bgzipWriter.WriteLine(header);
                }
                var writer = new BedGraphWriter(bgzipWriter);
                writer.WriteLines(bedGraphEntries);
            }
        }
    }
}