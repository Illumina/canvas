using System.Collections.Generic;
using System.Runtime.InteropServices.ComTypes;
using CanvasCommon;
using Xunit;

namespace CanvasTest
{
    public class CanvasSegmentTests
    {
        CanvasSegment _segment = new CanvasSegment("ada", 1000, 2000, new List<SampleGenomicBin>());

        [Fact]
        public void GetAltCopyNumbers_cnvType_is_reference()
        {
            var altString = _segment.GetAltCopyNumbers(CnvType.Reference);
            Assert.Equal(".", altString);
        }

        [Fact]
        public void GetAltCopyNumbers_major_chrom_count_null()
        {
            _segment.MajorChromosomeCount = null;
            var altString1 = _segment.GetAltCopyNumbers(CnvType.ComplexCnv);
            var altString2 = _segment.GetAltCopyNumbers(CnvType.Gain);
            var altString3 = _segment.GetAltCopyNumbers(CnvType.Loss);
            var altString4 = _segment.GetAltCopyNumbers(CnvType.LossOfHeterozygosity);
            Assert.Equal("<CNV>", altString1);
            Assert.Equal("<CNV>", altString2);
            Assert.Equal("<CNV>", altString3);
            Assert.Equal("<CNV>", altString4);
        }

        [Fact]
        public void GetAltCopyNumbers_one_allele_is_ref()
        {
            _segment.MajorChromosomeCount = 3;
            _segment.CopyNumber = 4;
            var altString1 = _segment.GetAltCopyNumbers(CnvType.ComplexCnv);
            var altString2 = _segment.GetAltCopyNumbers(CnvType.Gain);
            var altString3 = _segment.GetAltCopyNumbers(CnvType.Loss);
            var altString4 = _segment.GetAltCopyNumbers(CnvType.LossOfHeterozygosity);
            Assert.Equal("<CN3>", altString1);
            Assert.Equal("<CN3>", altString2);
            Assert.Equal("<CN3>", altString3);
            Assert.Equal("<CN3>", altString4);
        }

        [Fact]
        public void GetAltCopyNumbers_none_allele_is_ref()
        {
            _segment.MajorChromosomeCount = 3;
            _segment.CopyNumber = 5;
            var altString1 = _segment.GetAltCopyNumbers(CnvType.ComplexCnv);
            var altString2 = _segment.GetAltCopyNumbers(CnvType.Gain);
            var altString3 = _segment.GetAltCopyNumbers(CnvType.Loss);
            var altString4 = _segment.GetAltCopyNumbers(CnvType.LossOfHeterozygosity);
            Assert.Equal("<CN2>,<CN3>", altString1);
            Assert.Equal("<CN2>,<CN3>", altString2);
            Assert.Equal("<CN2>,<CN3>", altString3);
            Assert.Equal("<CN2>,<CN3>", altString4);
        }
    }
}
