using System;
using System.Collections.Generic;
using CanvasCommon;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace CanvasTest
{
    [TestClass]
    public class TestSegments
    {


        [TestMethod]
        public void TestCIPOS()
        {
            // Merge two segments, and confirm we keep the correct confidence intervals post-merge:
            List<float> counts = new List<float>() { 100, 90, 110, 100, 95, 105 };
            CanvasSegment segment = new CanvasSegment("chr1", 1245, 678910, counts);
            segment.StartConfidenceInterval = new Tuple<int, int>(-100, 100);
            segment.EndConfidenceInterval = new Tuple<int, int>(-80, 80);
            CanvasSegment segment2 = new CanvasSegment("chr1", 678910, 8787888, counts);
            segment2.StartConfidenceInterval = new Tuple<int, int>(-50, 50);
            segment2.EndConfidenceInterval = new Tuple<int, int>(-30, 30);
            segment.MergeIn(segment2);
            Assert.AreEqual(segment.End, 8787888);
            Assert.AreEqual(segment.EndConfidenceInterval.Item1, -30);
            Assert.AreEqual(segment.StartConfidenceInterval.Item2, 100);
        }

        [TestMethod]
        public void TestBins()
        {
            GenomicBin bin = new GenomicBin("chr1", 12345, 678910, 20, 100);
            Assert.AreEqual(bin.Size, 666565);
        }

        [TestMethod]
        public void TestSegment()
        {
            List<float> counts = new List<float>() {100, 90, 110, 100, 95, 105};
            CanvasSegment seg1 = new CanvasSegment("chr17", 100000000, 110000000, counts);
            // Silly constructor tests:
            Assert.AreEqual(seg1.Begin, 100000000);
            Assert.AreEqual(seg1.End, 110000000);
            Assert.AreEqual(seg1.BinCount, counts.Count);
            Assert.AreEqual(seg1.Chr, "chr17");
            // Property test:
            Assert.AreEqual(seg1.MeanCount, 100, 0.01);

            // Build a second segment, and merge them, and test results:
            CanvasSegment seg2 = new CanvasSegment("chr17", 110000000, 120000000, counts);
            seg1.MergeIn(seg2);
            Assert.AreEqual(seg1.Counts.Count, 12);
            Assert.AreEqual(seg1.End, seg2.End);
        }

        [TestMethod]
        public void TestSegmentStats()
        {
            List<float> counts = new List<float>() { 80, 79, 78, 77, 2 };
            List<CanvasSegment> segments = new List<CanvasSegment>();
            for (int index = 0; index < 10; index++)
            {
                CanvasSegment seg = new CanvasSegment("chr10", 1000000 * index, 1000000 * (index + 1), counts);
                segments.Add(seg);
            }
            double expectedCount = CanvasSegment.ExpectedCount(segments);
            Assert.AreEqual(expectedCount, 78, 0.01);
        }

        [TestMethod]
        public void TestMergeSegments()
        {
            // Construct several segments, and invoke CanvasSegment.MergeSegments, and ensure that the expected
            // merges (and no others) occurred.
            List<CanvasSegment> allSegments = new List<CanvasSegment>();
            List<float> counts = new List<float>();
            // Chr1 gets five segments and we should merge to three:
            CanvasSegment seg = new CanvasSegment("chr1", 1000000, 2000000, counts);
            seg.CopyNumber = 2;
            allSegments.Add(seg);
            seg = new CanvasSegment("chr1", 2000000, 2000100, counts);
            seg.CopyNumber = 3;
            allSegments.Add(seg);
            seg = new CanvasSegment("chr1", 2000100, 3000000, counts);
            seg.CopyNumber = 2;
            allSegments.Add(seg);
            seg = new CanvasSegment("chr1", 3000000, 3100000, counts);
            seg.CopyNumber = 3;
            allSegments.Add(seg);
            seg = new CanvasSegment("chr1", 3100000, 4000000, counts);
            seg.CopyNumber = 2;
            allSegments.Add(seg);

            // Chr2 gets segments with a large gap between, so can't merge:
            seg = new CanvasSegment("chr2", 1000000, 2000000, counts);
            seg.CopyNumber = 2;
            allSegments.Add(seg);
            seg = new CanvasSegment("chr2", 3000000, 3000100, counts);
            seg.CopyNumber = 3;
            allSegments.Add(seg);
            seg = new CanvasSegment("chr2", 4000000, 5000000, counts);
            seg.CopyNumber = 2;
            allSegments.Add(seg);

            // Chr3 has three segments that all merge to 1 big one:
            seg = new CanvasSegment("chr3", 1000000, 2000000, counts);
            seg.CopyNumber = 2;
            allSegments.Add(seg);
            seg = new CanvasSegment("chr3", 2000000, 3000000, counts);
            seg.CopyNumber = 2;
            allSegments.Add(seg);
            seg = new CanvasSegment("chr3", 3000000, 4000000, counts);
            seg.CopyNumber = 2;
            allSegments.Add(seg);

            CanvasSegment.MergeSegments(ref allSegments, 50000, 10000);
            Dictionary<string, List<CanvasSegment>> segmentsByChromosome = CanvasSegment.GetSegmentsByChromosome(allSegments);
            Assert.AreEqual(segmentsByChromosome["chr1"].Count, 3);
            Assert.AreEqual(segmentsByChromosome["chr2"].Count, 3);
            Assert.AreEqual(segmentsByChromosome["chr3"].Count, 1);
        }

    }
}
