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

        [Fact]
        public void GetAltAllelesAndGenotypes_various_genotypes()
        {
            int[][] alleleCopyNumbers =
            {
                new[] { -1, int.MaxValue}, // <DUP>
                new[] {1, 1 }, // Ref
                new[] {0, 1 }, // Loss
                new[] {1,2}, // Gain
                new[] {0,3}, // Gain
                new[] {2,2}, // Gain
                new [] {0, 2}, // LOH
                new [] {-1, -1}, // CN=2, MCC null
                new [] {-1 } // Ploidy == 0
            };
            var (altAlleleString, sampleGenotypes) = CanvasSegmentWriter.GetAltAllelesAndGenotypes(alleleCopyNumbers);

            string[] expectedGenotypes = { "./4", "0/0", "0/1", "0/2", "1/3", "2/2","1/2", "./.", "." };
            Assert.Equal("<CN0>,<CN2>,<CN3>,<DUP>", altAlleleString);
            Assert.Equal(expectedGenotypes, sampleGenotypes);
        }

        [Fact]
        public void GetAltAllelesAndGenotypes_only_reference()
        {
            int[][] alleleCopyNumbers =
            {
                new[] {1, 1 }, // Ref
                new[] {1, 1 } // Ref again
            };
            var (altAlleleString, sampleGenotypes) = CanvasSegmentWriter.GetAltAllelesAndGenotypes(alleleCopyNumbers);

            string[] expectedGenotypes = { "0/0", "0/0" };
            Assert.Equal(".", altAlleleString);
            Assert.Equal(expectedGenotypes, sampleGenotypes);
        }

        [Fact]
        public void GetAltAllelesAndGenotypes_hemizygous_regions()
        {
            int[][] alleleCopyNumbers =
            {
                new[] {0}, // Loss
                new[] {1}, // Ref
                new[] {2} // Gain
            };

            var (altAlleleString, sampleGenotypes) = CanvasSegmentWriter.GetAltAllelesAndGenotypes(alleleCopyNumbers);

            string[] expectedGenotypes = { "1", "0", "2" };
            Assert.Equal("<CN0>,<CN2>", altAlleleString);
            Assert.Equal(expectedGenotypes, sampleGenotypes);
        }

        [Fact]
        public void GetAltAllelesAndGenotypes_one_sample_ploidy_zero()
        {
            int[][] alleleCopyNumbers =
            {
                new[] {-1 }
            };
            var (altAlleleString, sampleGenotypes) = CanvasSegmentWriter.GetAltAllelesAndGenotypes(alleleCopyNumbers);

            string[] expectedGenotypes = { "." };
            Assert.Equal(".", altAlleleString);
            Assert.Equal(expectedGenotypes, sampleGenotypes);
        }
    }
}
