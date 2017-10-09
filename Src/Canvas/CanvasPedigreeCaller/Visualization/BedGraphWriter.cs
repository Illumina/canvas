using System.Collections.Generic;
using Illumina.Common;
using Isas.SequencingFiles;
using Isas.SequencingFiles.Bed;

namespace CanvasPedigreeCaller.Visualization
{
    /// <summary> 
    /// A class to write out lines in bedgraph format
    /// Decimal values are written out with at most 29 digits of precision (trailing zeros ommitted)
    /// Any necessary rounding should be done on the BedGraphEntry objects before writing them out
    /// </summary>
    public class BedGraphWriter
    {
        private readonly BgzipOrStreamWriter _writer;

        public BedGraphWriter(BgzipOrStreamWriter writer)
        {
            _writer = writer;
        }

        public void WriteLine(BedGraphEntry entry)
        {
            var bedGraphLine = GetBedGraphLine(entry);
            _writer.WriteLine(bedGraphLine);
        }

        private string GetBedGraphLine(BedGraphEntry entry)
        {
            return $"{entry.Chromosome}\t{entry.Interval.Start}\t{entry.Interval.End}\t{GetValueAsString(entry.Value)}";
        }

        private string GetValueAsString(decimal value)
        {
            return value.ToString($"0.{new string('#', 29)}");
        }
    }
    
    public static class BedGraphWriterExtensions
    {
        public static void WriteLines(this BedGraphWriter writer, IEnumerable<BedGraphEntry> entries)
        {
            entries.ForEach(writer.WriteLine);
        }
    }
}