using System.Collections.Generic;
using CanvasCommon;
using Illumina.Common.FileSystem;
using Isas.SequencingFiles;
using Isas.SequencingFiles.Bed;

namespace CanvasPedigreeCaller.Visualization
{
    public interface ICoverageBedGraphWriter
    {
        void Write(IReadOnlyList<CanvasSegment> segments, IFileLocation location);
    }

    public class NormalizedCoverageWriterFacade : ICoverageBedGraphWriter
    {
        private readonly NormalizedCoverageBedGraphWriter _normalizedCoverageBedGraphWriter;

        public NormalizedCoverageWriterFacade(NormalizedCoverageBedGraphWriter normalizedCoverageCalculator)
        {
            _normalizedCoverageBedGraphWriter = normalizedCoverageCalculator;
        }

        public void Write(IReadOnlyList<CanvasSegment> segments, IFileLocation location)
        {
            using (var bgzipWriter = new BgzipOrStreamWriter(location.FullName))
            {
                var writer = new BedGraphWriter(bgzipWriter);
                _normalizedCoverageBedGraphWriter.Write(segments, writer);
            }
        }
    }
}