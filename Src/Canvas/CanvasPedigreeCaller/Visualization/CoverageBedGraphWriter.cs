using System.Collections.Generic;
using CanvasCommon;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Isas.ClassicBioinfoTools.KentUtils;

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