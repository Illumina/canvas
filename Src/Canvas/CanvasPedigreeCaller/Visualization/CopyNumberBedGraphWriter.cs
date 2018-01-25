using System.Collections.Generic;
using CanvasCommon;
using Illumina.Common.FileSystem;
using Isas.SequencingFiles.Bed;

namespace CanvasPedigreeCaller.Visualization
{
    public class CopyNumberBedGraphWriter : ICopyNumberBedGraphWriter
    {
        private readonly CopyNumberBedGraphCalculator _calculator;

        public CopyNumberBedGraphWriter(CopyNumberBedGraphCalculator calculator)
        {
            _calculator = calculator;
        }

        public void Write(IReadOnlyList<CanvasSegment> segments, PloidyInfo ploidyInfo, IFileLocation location)
        {
            _calculator.Calculate(segments, ploidyInfo);
            var calculator = new CopyNumberCalculatorAdapter(_calculator, ploidyInfo);
            var writer = new BedGraphWriterFacade(calculator);
            writer.Write(segments, location);
        }

        private class CopyNumberCalculatorAdapter : IBedGraphCalculator
        {
            private readonly CopyNumberBedGraphCalculator _calculator;
            private readonly PloidyInfo _ploidyInfo;

            public CopyNumberCalculatorAdapter(CopyNumberBedGraphCalculator calculator, PloidyInfo ploidyInfo)
            {
                _calculator = calculator;
                _ploidyInfo = ploidyInfo;
            }

            public IEnumerable<BedGraphEntry> Calculate(IReadOnlyList<CanvasSegment> segments)
            {
                return _calculator.Calculate(segments, _ploidyInfo);
            }
        }
    }
}