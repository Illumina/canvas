using CanvasCommon;
using CanvasCommon.Visualization;
using Illumina.Common.FileSystem;
using System.Collections.Generic;

namespace CanvasPedigreeCaller.Visualization
{
    public class CoverageBedGraphWriter
    {
        private readonly IBedGraphWriter _writer;
        private readonly BaseNormalizedCoverageCalculator _bedGraphCalculator;

        public CoverageBedGraphWriter(IBedGraphWriter writer, BaseNormalizedCoverageCalculator bedGraphCalculator)
        {
            _writer = writer;
            _bedGraphCalculator = bedGraphCalculator;
        }

        public void Write(IReadOnlyList<CanvasSegment> segments, IFileLocation location, double normalizationFactor, string header = null)
        {
            var entries = _bedGraphCalculator.Calculate(segments, normalizationFactor);
            _writer.Write(entries, location, header);
        }
    }
}