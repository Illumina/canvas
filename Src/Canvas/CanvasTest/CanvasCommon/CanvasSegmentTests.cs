using System.Collections.Generic;
using CanvasCommon;
using Xunit;

namespace CanvasTest
{
    public class CanvasSegmentTests
    {
        private readonly CanvasSegment _segmentCn0 = new CanvasSegment("ada", 1000, 2000, new List<SampleGenomicBin>()) { CopyNumber = 0 };
        private readonly CanvasSegment _segmentCn1 = new CanvasSegment("ada", 1000, 2000, new List<SampleGenomicBin>()) { CopyNumber = 1 };
        private readonly CanvasSegment _segmentCn2 = new CanvasSegment("ada", 1000, 2000, new List<SampleGenomicBin>()) { CopyNumber = 2 };
        private readonly CanvasSegment _segmentCn2Mcc2 = new CanvasSegment("ada", 1000, 2000, new List<SampleGenomicBin>()) { CopyNumber = 2, MajorChromosomeCount = 2};
        private readonly CanvasSegment _segmentCn3 = new CanvasSegment("ada", 1000, 2000, new List<SampleGenomicBin>()) {CopyNumber = 3} ;
        private readonly CanvasSegment _segmentCn3Mcc2 = new CanvasSegment("ada", 1000, 2000, new List<SampleGenomicBin>()) { CopyNumber = 3, MajorChromosomeCount = 2 };

        [Fact]
        public void GetCnvTypeAndAlleleCopyNumbers_referenceCn_is_two()
        {
            int referenceCopyNumber = 2;
            var (cnvTypeCn0, alleleCopyNumbersCn0) = _segmentCn0.GetCnvTypeAndAlleleCopyNumbers(referenceCopyNumber);
            var (cnvTypeCn1, alleleCopyNumbersCn1) = _segmentCn1.GetCnvTypeAndAlleleCopyNumbers(referenceCopyNumber);
            var (cnvTypeCn2, alleleCopyNumbersCn2) = _segmentCn2.GetCnvTypeAndAlleleCopyNumbers(referenceCopyNumber);
            var (cnvTypeCn2Mcc2, alleleCopyNumbersCn2Mcc2) = _segmentCn2Mcc2.GetCnvTypeAndAlleleCopyNumbers(referenceCopyNumber);
            var (cnvTypeCn3, alleleCopyNumbersCn3) = _segmentCn3.GetCnvTypeAndAlleleCopyNumbers(referenceCopyNumber);
            var (cnvTypeCn3Mcc2, alleleCopyNumbersCn3Mcc2) = _segmentCn3Mcc2.GetCnvTypeAndAlleleCopyNumbers(referenceCopyNumber);

            Assert.Equal(CnvType.Loss, cnvTypeCn0);
            Assert.Equal(new [] {0, 0}, alleleCopyNumbersCn0);
            Assert.Equal(CnvType.Loss, cnvTypeCn1);
            Assert.Equal(new[] { 0, 1 }, alleleCopyNumbersCn1);
            Assert.Equal(CnvType.Reference, cnvTypeCn2);
            Assert.Equal(new[] { -1, -1 }, alleleCopyNumbersCn2);
            Assert.Equal(CnvType.LossOfHeterozygosity, cnvTypeCn2Mcc2);
            Assert.Equal(new[] { 0, 2 }, alleleCopyNumbersCn2Mcc2);
            Assert.Equal(CnvType.Gain, cnvTypeCn3);
            Assert.Equal(new[] { -1, int.MaxValue}, alleleCopyNumbersCn3);
            Assert.Equal(CnvType.Gain, cnvTypeCn3Mcc2);
            Assert.Equal(new[] { 1, 2 }, alleleCopyNumbersCn3Mcc2);
        }

        [Fact]
        public void GetCnvTypeAndAlleleCopyNumbers_referenceCn_is_one()
        {
            int referenceCopyNumber = 1;
            var (cnvTypeCn0, alleleCopyNumbersCn0) = _segmentCn0.GetCnvTypeAndAlleleCopyNumbers(referenceCopyNumber);
            var (cnvTypeCn1, alleleCopyNumbersCn1) = _segmentCn1.GetCnvTypeAndAlleleCopyNumbers(referenceCopyNumber);
            var (cnvTypeCn2, alleleCopyNumbersCn2) = _segmentCn2.GetCnvTypeAndAlleleCopyNumbers(referenceCopyNumber);

            Assert.Equal(CnvType.Loss, cnvTypeCn0);
            Assert.Equal(new[] { 0 }, alleleCopyNumbersCn0);
            Assert.Equal(CnvType.Reference, cnvTypeCn1);
            Assert.Equal(new[] { 1 }, alleleCopyNumbersCn1);
            Assert.Equal(CnvType.Gain, cnvTypeCn2);
            Assert.Equal(new[] { 2 }, alleleCopyNumbersCn2);
        }

        [Fact]
        public void GetCnvTypeAndAlleleCopyNumbers_referenceCn_is_zero()
        {
            int referenceCopyNumber = 0;
            var (cnvTypeCn0, alleleCopyNumbersCn0) = _segmentCn0.GetCnvTypeAndAlleleCopyNumbers(referenceCopyNumber);
            var (cnvTypeCn1, alleleCopyNumbersCn1) = _segmentCn1.GetCnvTypeAndAlleleCopyNumbers(referenceCopyNumber);
            var (cnvTypeCn2, alleleCopyNumbersCn2) = _segmentCn2.GetCnvTypeAndAlleleCopyNumbers(referenceCopyNumber);

            Assert.Equal(CnvType.Reference, cnvTypeCn0);
            Assert.Equal(new[] { -1 }, alleleCopyNumbersCn0);
            Assert.Equal(CnvType.Gain, cnvTypeCn1);
            Assert.Equal(new[] { -1 }, alleleCopyNumbersCn1);
            Assert.Equal(CnvType.Gain, cnvTypeCn2);
            Assert.Equal(new[] { -1 }, alleleCopyNumbersCn2);
        }
    }
}
