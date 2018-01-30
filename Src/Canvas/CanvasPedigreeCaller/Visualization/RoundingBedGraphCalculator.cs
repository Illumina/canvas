using System;
using System.Collections.Generic;
using System.Linq;
using CanvasCommon;
using Isas.SequencingFiles.Bed;

namespace CanvasPedigreeCaller.Visualization
{

    public class RoundingBedGraphCalculator : IBedGraphCalculator
    {
        private readonly IBedGraphCalculator _calculator;
        private readonly int _fractionalDigits;

        public RoundingBedGraphCalculator(IBedGraphCalculator calculator, int fractionalDigits)
        {
            _calculator = calculator;
            _fractionalDigits = fractionalDigits;
        }

        public IEnumerable<BedGraphEntry> Calculate(IReadOnlyList<CanvasSegment> segments)
        {
            var entries = _calculator.Calculate(segments);
            return entries.Select(GetRoundedBedGraphEntry);
        }

        private BedGraphEntry GetRoundedBedGraphEntry(BedGraphEntry entry)
        {
            var roundedCoverage = Math.Round(entry.Value, _fractionalDigits, MidpointRounding.AwayFromZero);
            return new BedGraphEntry(entry.Chromosome, entry.Interval, roundedCoverage);
        }
    }
}