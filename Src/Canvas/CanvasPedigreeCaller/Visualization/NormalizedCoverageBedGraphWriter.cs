using System;
using System.Collections.Generic;
using System.Linq;
using CanvasCommon;
using Isas.SequencingFiles.Bed;

namespace CanvasPedigreeCaller.Visualization
{
    public class NormalizedCoverageBedGraphWriter
    {
        private readonly NormalizedCoverageCalculator _normalizedCoverageCalculator;
        private readonly int _fractionalDigits;

        public NormalizedCoverageBedGraphWriter(NormalizedCoverageCalculator normalizedCoverageCalculator, int fractionalDigits)
        {
            _normalizedCoverageCalculator = normalizedCoverageCalculator;
            _fractionalDigits = fractionalDigits;
        }
        public void Write(IReadOnlyList<CanvasSegment> segments, BedGraphWriter bedGraphWriter)
        {
            var entries = _normalizedCoverageCalculator.Calculate(segments);
            entries = entries.Select(GetRoundedBedGraphEntry);
            bedGraphWriter.WriteLines(entries);
        }

        private BedGraphEntry GetRoundedBedGraphEntry(BedGraphEntry entry)
        {
            var roundedCoverage = Math.Round(entry.Value, _fractionalDigits, MidpointRounding.AwayFromZero);
            return new BedGraphEntry(entry.Chromosome, entry.Interval, roundedCoverage);
        }
    }
}