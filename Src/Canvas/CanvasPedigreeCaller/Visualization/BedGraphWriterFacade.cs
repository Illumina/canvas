using System;
using System.Collections.Generic;
using CanvasCommon;
using Illumina.Common.FileSystem;
using Isas.SequencingFiles;
using Isas.SequencingFiles.Bed;

namespace CanvasPedigreeCaller.Visualization
{
    public class BedGraphWriterFacade : IBedGraphWriter
    {
        private readonly IBedGraphCalculator _bedGraphCalculator;

        public BedGraphWriterFacade(IBedGraphCalculator bedGraphCalculator)
        {
            _bedGraphCalculator = bedGraphCalculator;
        }

        public void Write(IReadOnlyList<CanvasSegment> segments, IFileLocation location)
        {
            using (var bgzipWriter = new BgzipOrStreamWriter(location.FullName))
            {
                var writer = new BedGraphWriter(bgzipWriter);
                var entries = _bedGraphCalculator.Calculate(segments);
                writer.WriteLines(entries);
            }
        }
    }
}