using System.Collections.Generic;
using System.IO;
using System.Linq;
using CanvasCommon.Visualization;
using Illumina.Common.CSV;
using Illumina.Common.FileSystem;
using Isas.Framework.DataTypes;
using Isas.SequencingFiles;
using Isas.SequencingFiles.Bed;

namespace Canvas.Visualization
{
    public class BAlleleBedGraphWriter : IBAlleleBedGraphWriter
    {
        private readonly IBgzfBedGraphWriter _bedGraphWriter;

        public BAlleleBedGraphWriter(IBgzfBedGraphWriter bedGraphWriter)
        {
            _bedGraphWriter = bedGraphWriter;
        }

        public void Write(IFileLocation bafFile, BgzfFile bAllelesFile)
        {
            var bedGraphEntries = GetAlleleFrequencies(bafFile);
            _bedGraphWriter.Write(bedGraphEntries, bAllelesFile);
        }

        private IEnumerable<BedGraphEntry> GetAlleleFrequencies(IFileLocation bafFile)
        {
            return File.ReadLines(bafFile.FullName).Select(GetAlleleFrequencyEntry);
        }

        private BedGraphEntry GetAlleleFrequencyEntry(string bafLine)
        {
            var bafFields = CSVReader.ParseCommaDelimitedLine(bafLine);
            var chromosome = bafFields[0];
            var oneBasedPosition = int.Parse(bafFields[1]);
            var alleleFrequency = decimal.Parse(bafFields[2]);
            var bedPosition = new BedInterval(oneBasedPosition - 1, oneBasedPosition);
            return new BedGraphEntry(chromosome, bedPosition, alleleFrequency);
        }
    }
}