﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using CanvasCommon;
using Isas.Framework.Logging;
using Isas.SequencingFiles;
using Isas.SequencingFiles.Vcf;
using JetBrains.Annotations;
using Xunit;

namespace CanvasTest
{
    public class TestSegments
    {


        [Fact]
        public void TestCIPOS()
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
            CanvasSegment segment = new CanvasSegment("chr1", 1245, 678910, counts);
            segment.StartConfidenceInterval = new Tuple<int, int>(-100, 100);
            segment.EndConfidenceInterval = new Tuple<int, int>(-80, 80);
            CanvasSegment segment2 = new CanvasSegment("chr1", 678910, 8787888, counts);
            segment2.StartConfidenceInterval = new Tuple<int, int>(-50, 50);
            segment2.EndConfidenceInterval = new Tuple<int, int>(-30, 30);
            segment.MergeIn(segment2);
            Assert.Equal(segment.End, 8787888);
            Assert.Equal(segment.EndConfidenceInterval.Item1, -30);
            Assert.Equal(segment.StartConfidenceInterval.Item2, 100);
        }

        [Fact]
        public void TestBins()
        {
            SampleGenomicBin bin = new SampleGenomicBin("chr1", 12345, 678910, 20, 100);
            Assert.Equal(bin.Size, 666565);
        }

        [Fact]
        public void TestSegment()
        {
            var counts = new List<SampleGenomicBin>
            {
                new SampleGenomicBin("chr17", 100000000, 110000000, 0, 100),
                new SampleGenomicBin("chr17", 100000000, 110000000, 0, 90),
                new SampleGenomicBin("chr17", 100000000, 110000000, 0, 110),
                new SampleGenomicBin("chr17", 100000000, 110000000, 0, 100),
                new SampleGenomicBin("chr17", 100000000, 110000000, 0, 95),
                new SampleGenomicBin("chr17", 100000000, 110000000, 0, 105)
            };
            var seg1 = new CanvasSegment("chr17", 100000000, 110000000, counts);
            // Silly constructor tests:
            Assert.Equal(seg1.Begin, 100000000);
            Assert.Equal(seg1.End, 110000000);
            Assert.Equal(seg1.BinCount, counts.Count);
            Assert.Equal(seg1.Chr, "chr17");
            // Property test:
            Assert.Equal(seg1.MeanCount, 100, 2);

            // Build a second segment, and merge them, and test results:
            var seg2 = new CanvasSegment("chr17", 110000000, 120000000, counts);
            seg1.MergeIn(seg2);
            Assert.Equal(seg1.Counts.Count, 12);
            Assert.Equal(seg1.End, seg2.End);
        }

        [Fact]
        public void TestSegmentStats()
        {
            var counts = new List<SampleGenomicBin>
            {
                new SampleGenomicBin("chr10", 1000000, 1000001, 0, 80),
                new SampleGenomicBin("chr10", 1000000, 1000001, 0, 79),
                new SampleGenomicBin("chr10", 1000000, 1000001, 0, 78),
                new SampleGenomicBin("chr10", 1000000, 1000001, 0, 77),
                new SampleGenomicBin("chr10", 1000000, 1000001, 0, 2)
            };
            var segments = new List<CanvasSegment>();
            for (int index = 0; index < 10; index++)
            {
                var seg = new CanvasSegment("chr10", 1000000 * index, 1000000 * (index + 1), counts);
                segments.Add(seg);
            }
            double expectedCount = CanvasSegment.ExpectedCount(segments);
            Assert.Equal(expectedCount, 78, 2);
        }

        [Fact]
        public void TestMergeSegments()
        {
            // Construct several segments, and invoke CanvasSegment.MergeSegments, and ensure that the expected
            // merges (and no others) occurred.
            List<CanvasSegment> allSegments = new List<CanvasSegment>();
            List<SampleGenomicBin> counts = new List<SampleGenomicBin>();
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
            var segmentsByChromosome = CanvasSegment.GetSegmentsByChromosome(allSegments);
            Assert.Equal(segmentsByChromosome["chr1"].Count, 3);
            Assert.Equal(segmentsByChromosome["chr2"].Count, 3);
            Assert.Equal(segmentsByChromosome["chr3"].Count, 1);
        }

        [Fact]
        public void TestReadSegments()
        {
            var partitioned = "";
            partitioned += "chr22\t1\t10\t14.00\t0\n";
            partitioned += "chr22\t10\t30\t31.00\t1\n";
            partitioned += "chr22\t30\t40\t6.00\t2\n";
            var stringReader = new StringReader(partitioned);
            Segments segments;
            using (var reader = new GzipOrTextReader(stringReader))
            {
                segments = Segments.ReadSegments(reader);
            }

            Assert.Equal(segments.GetSegmentsForChromosome("chr22"), segments.AllSegments);
            var confidenceInterval = segments.AllSegments.First().StartConfidenceInterval;
            AssertConfidenceInterval(-5, 5, confidenceInterval);

            confidenceInterval = segments.AllSegments.First().EndConfidenceInterval;
            AssertConfidenceInterval(-5, 10, confidenceInterval);

            confidenceInterval = segments.AllSegments[1].StartConfidenceInterval;
            AssertConfidenceInterval(-5, 10, confidenceInterval);

            confidenceInterval = segments.AllSegments[1].EndConfidenceInterval;
            AssertConfidenceInterval(-10, 5, confidenceInterval);

            confidenceInterval = segments.AllSegments.Last().StartConfidenceInterval;
            AssertConfidenceInterval(-10, 5, confidenceInterval);

            confidenceInterval = segments.AllSegments.Last().EndConfidenceInterval;
            AssertConfidenceInterval(-5, 5, confidenceInterval);
        }

        [AssertionMethod]
        private void AssertConfidenceInterval(int expectedLower, int expectedUpper, Tuple<int, int> interval)
        {
            Assert.Equal(expectedLower, interval.Item1);
            Assert.Equal(expectedUpper, interval.Item2);
        }

        [Fact]
        public void TestReadFrequencies()
        {
            var intervals = new List<BedInterval>
            {
                new BedInterval(1, 50),
                new BedInterval(51, 150),
            };
            const string chr = "chr22";
            var intervalsByChromosome = new Dictionary<string, List<BedInterval>> {{chr, intervals}};
            var variantCounts = "";
            variantCounts += "chr22\t10\tC\tT\t20\t10\n";
            variantCounts += "chr22\t20\tC\tT\t30\t20\n";
            variantCounts += "chr22\t100\tC\tT\t40\t30\n";
            var stringReader = new StringReader(variantCounts);
            using (var reader = new GzipOrTextReader(stringReader))
            {
                Dictionary<string, List<Balleles>>  allelesByChromosome = 
                    CanvasIO.ReadFrequencies(reader, intervalsByChromosome);
                Assert.Equal(allelesByChromosome[chr].Count, intervals.Count);
                Assert.Equal(2, allelesByChromosome[chr].First().Size());
                Assert.Equal(1, allelesByChromosome[chr].Last().Size());
            }
        }

    }
}
