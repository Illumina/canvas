using Xunit;
using CanvasSNV;
using Isas.SequencingFiles.Vcf;

namespace CanvasTest
{
    public class TestCanvasSNV
    {
        [Theory]
        [InlineData("A", "T", 0, 0, null)]
        [InlineData("A", "T", 1, 3, 0.25)]
        [InlineData("T", "A", 1, 3, 0.75)]
        [InlineData("T", "G", 1, 3, 0.25)]
        [InlineData("G", "C", 1, 3, 0.25)]
        [InlineData("A", "A", 1, 3, 0.75)]
        public void TestGetBAlleleFrequency(string refAllele, string altAllele, int refCount, int altCount,
            double? expectedFreq)
        {
            VcfVariant variant = new VcfVariant();
            variant.ReferenceAllele = refAllele;
            variant.VariantAlleles = new string[] { altAllele };

            double? freq = SNVReviewer.GetBAlleleFrequency(variant, refCount, altCount);
            Assert.Equal(expectedFreq, freq);
        }
    }
}
