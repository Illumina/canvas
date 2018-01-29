using System.Collections.Generic;
using CanvasCommon;
using CanvasCommon.Visualization;
using Isas.Framework.DataTypes;

namespace CanvasPedigreeCaller.Visualization
{
    public class CopyNumberBedGraphWriter : ICopyNumberBedGraphWriter
    {
        private readonly IBgzfBedGraphWriter _writer;
        private readonly CopyNumberBedGraphCalculator _calculator;

        public CopyNumberBedGraphWriter(IBgzfBedGraphWriter writer, CopyNumberBedGraphCalculator calculator)
        {
            _writer = writer;
            _calculator = calculator;
        }

        public void Write(IReadOnlyList<CanvasSegment> segments, PloidyInfo ploidyInfo, BgzfFile location)
        {
            var entries = _calculator.Calculate(segments, ploidyInfo);
            _writer.Write(entries, location);
        }
    }
}