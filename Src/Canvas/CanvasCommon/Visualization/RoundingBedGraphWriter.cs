using System;
using System.Collections.Generic;
using System.Linq;
using Illumina.Common.FileSystem;
using Isas.SequencingFiles.Bed;

namespace CanvasCommon.Visualization
{
   
    public class RoundingBedGraphWriter : IBedGraphWriter
    {
        private readonly IBedGraphWriter _bedGraphWriter;
        private readonly int _fractionalDigits;

        public RoundingBedGraphWriter(IBedGraphWriter bedGraphWriter, int fractionalDigits)
        {
            _bedGraphWriter = bedGraphWriter;
            _fractionalDigits = fractionalDigits;
        }

        public void Write(IEnumerable<BedGraphEntry> bedGraphEntries, IFileLocation location, string header=null)
        {
            var roundedEntries = bedGraphEntries.Select(GetRoundedBedGraphEntry);
            _bedGraphWriter.Write(roundedEntries, location, header);
        }

        private BedGraphEntry GetRoundedBedGraphEntry(BedGraphEntry entry)
        {
            var roundedCoverage = Math.Round(entry.Value, _fractionalDigits, MidpointRounding.AwayFromZero);
            return new BedGraphEntry(entry.Chromosome, entry.Interval, roundedCoverage);
        }
    }
}