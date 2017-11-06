using System.IO;
using CanvasCommon;
using Isas.SequencingFiles;
using Xunit;

namespace CanvasTest
{
    public class CanvasSegmentWriterTests
    {
        [Fact]
        public void WriteHeaderAllAltCnTags_output_is_expected()
        {
            const int maxCopyNum = 3;
            const string newLine = "\n";
            string expected = string.Join(newLine, "##ALT=<ID=CN0,Description=\"Copy number allele: 0 copies\">", "##ALT=<ID=CN2,Description=\"Copy number allele: 2 copies\">", "##ALT=<ID=CN3,Description=\"Copy number allele: 3 copies\">") + newLine;
            string output;
            using (var stringWriter = new StringWriter())
            {
                stringWriter.NewLine = newLine;
                CanvasSegmentWriter.WriteHeaderAllAltCnTags(new BgzipOrStreamWriter(stringWriter), maxCopyNum);
                output = stringWriter.ToString();
            }
            Assert.Equal(expected, output);
        }
    }
}
