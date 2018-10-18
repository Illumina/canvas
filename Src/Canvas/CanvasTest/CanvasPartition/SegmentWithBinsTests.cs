using System.Linq;
using CanvasPartition.Models;
using Xunit;

namespace CanvasTest.CanvasPartition
{
    public static class SegmentTestHelpers
    {
        public static void CheckSegment(SegmentWithBins segment, int expectedStart, int expectedEnd,
            double expectedMedCoverage, int expectedCount)
        {
            Assert.Equal((uint)expectedStart, segment.Start);
            Assert.Equal((uint)expectedEnd, segment.End);
            Assert.Equal(expectedMedCoverage, segment.MedianCoverage);
            Assert.Equal(expectedCount, segment.Bins.Count);

        }

    }
    public class SegmentWithBinsTests
    {
        [Fact]
        public void AddBinTest()
        {
            var bin1 = new Bin(100, 2000, 10d);
            var bin2 = new Bin(2500, 3000, 5d);
            var bin3 = new Bin(5000, 8000, 45d);

            var segmentWithBins = new SegmentWithBins(1, bin1);
            SegmentTestHelpers.CheckSegment(segmentWithBins, 100, 2000, 10, 1);

            segmentWithBins.AddBin(bin2);
            SegmentTestHelpers.CheckSegment(segmentWithBins, 100, 3000, 7.5, 2);

            segmentWithBins.AddBin(bin3);
            SegmentTestHelpers.CheckSegment(segmentWithBins, 100, 8000, 10, 3);

            // Order they are added in should not affect final result
            segmentWithBins = new SegmentWithBins(1, bin1);
            segmentWithBins.AddBin(bin3);
            segmentWithBins.AddBin(bin2);
            SegmentTestHelpers.CheckSegment(segmentWithBins, 100, 8000, 10, 3);

        }


    }
}