using CanvasCommon;
using CanvasCommon.Visualization;
using Illumina.Common.FileSystem;
using System.Collections.Generic;

namespace CanvasPedigreeCaller.Visualization
{
    public class CoverageBedGraphWriter
    {
        private readonly IBedGraphWriter _writer;
        private readonly NormalizedCoverageCalculator _bedGraphCalculator;

        public CoverageBedGraphWriter(IBedGraphWriter writer, NormalizedCoverageCalculator bedGraphCalculator)
        {
            _writer = writer;
            _bedGraphCalculator = bedGraphCalculator;
        }

        public void Write(IReadOnlyList<CanvasSegment> segments, IFileLocation location)
        {
            var entries = _bedGraphCalculator.Calculate(segments);
            _writer.Write(entries, location);
        }
    }
}