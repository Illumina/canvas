using System;
using System.Collections.Generic;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Isas.SequencingFiles;
using Isas.SequencingFiles.Bed;

namespace CanvasPedigreeCaller
{
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
            // decimal has at most 29 digits of precision
            // trailing zeros will not be output
            // any rounding should be done by the caller
            return value.ToString($"0.{new string('#', 29)}");
        }

        public static void WriteLines(IFileLocation location, IEnumerable<BedGraphEntry> entries)
        {
            using (var streamWriter = new BgzipOrStreamWriter(location.FullName))
            {
                var writer = new BedGraphWriter(streamWriter);
                writer.WriteLines(entries);
            }
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