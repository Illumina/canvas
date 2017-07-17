using CanvasCommon;
using System.Collections.Generic;
using Xunit;

namespace CanvasTest
{
    public class TestCommonCNVsSegments
    {
        [Fact]
        public void TestSplit()
        {
            // Merge two segments, and confirm we keep the correct confidence intervals post-merge:
            List<SampleGenomicBin> counts = new List<SampleGenomicBin>()
            {
                new SampleGenomicBin("chr1", 1, 2, 100),
                new SampleGenomicBin("chr1", 1, 2, 90),
                new SampleGenomicBin("chr1", 1, 2, 110),
                new SampleGenomicBin("chr1", 1, 2, 100),
                new SampleGenomicBin("chr1", 1, 2, 95),
                new SampleGenomicBin("chr1", 1, 2, 105)
            };
            var canvasSegments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 100000, 400000, counts),
            };
            var commonSegments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 200000, 300000, counts),
            };
            int index = 0;
            int defaultAlleleCountThreshold = 4;
            var haplotypeSegments = CanvasSegment.Split(canvasSegments, commonSegments, defaultAlleleCountThreshold, ref index, ref index);
            Assert.Equal(haplotypeSegments.SetA.Count, 1);
            Assert.Equal(haplotypeSegments.SetB.Count, 1);
        }
    }
}