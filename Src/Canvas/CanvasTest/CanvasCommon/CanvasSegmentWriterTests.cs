using System;
using System.Collections.Generic;
using System.IO;
using CanvasCommon;
using System.Text;
using Isas.SequencingFiles;
using Xunit;

namespace CanvasTest
{
    public class CanvasSegmentWriterTests
    {
        [Fact]
        public void WriteHeaderAllAltCnTags_output_is_expected()
        {
            int maxCopyNum = 3;
            string expected =
@"##ALT=<ID=CN0,Description=""Copy number allele: 0 copies"">
##ALT=<ID=CN2,Description=""Copy number allele: 2 copies"">
##ALT=<ID=CN3,Description=""Copy number allele: 3 copies"">
";
            string output;
            using (var stream = new MemoryStream())
            {
                using (var writer = new StreamWriter(stream))
                {
                    CanvasSegmentWriter.WriteHeaderAllAltCnTags(new BgzipOrStreamWriter(writer), maxCopyNum);
                }
                output = Encoding.UTF8.GetString(stream.ToArray());
            }
            Assert.Equal(expected, output);
        }
    }
}
