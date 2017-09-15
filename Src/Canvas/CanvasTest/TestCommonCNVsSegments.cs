using CanvasCommon;
using System.Collections.Generic;
using System.Linq;
using Xunit;

namespace CanvasTest
{
    public class TestCommonCNVsSegments
    {
        [Fact]
        public void TestSplitCommonCNVWithinCanvasCNV()
        {
            // scenario: common segment within Canvas segment
            // canvasSegment:   ----------------------------------
            // commonSegment:         -----------------

            var counts = new List<SampleGenomicBin>
            {
                new SampleGenomicBin("chr1", 100000, 100001, 100),
                new SampleGenomicBin("chr1", 150000, 150001, 90),
                new SampleGenomicBin("chr1", 200000, 200001, 110),
                new SampleGenomicBin("chr1", 250000, 250001, 100),
                new SampleGenomicBin("chr1", 300000, 300001, 95),
                new SampleGenomicBin("chr1", 350000, 350001, 105),
                new SampleGenomicBin("chr1", 400000, 400001, 105),
                new SampleGenomicBin("chr1", 450000, 450001, 105),
                new SampleGenomicBin("chr1", 500000, 500001, 105)
            };
            var canvasSegments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 100000, 500002, counts),
            };
            var commonSegments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 250000, 350001, counts.Skip(3).Take(3).ToList()),
            };
            var canvasSegmentsIndex = 0;
            var commonSegmentsIndex = 0;
            const int defaultReadCountsThreshold = 4;
            var haplotypeSegments = CanvasSegment.SplitCanvasSegments(canvasSegments, commonSegments, defaultReadCountsThreshold, ref canvasSegmentsIndex, ref commonSegmentsIndex);

            // transform into "haplotype" segments
            // canvasSegment:   ----------------------------------
            // commonSegment:   -----  ----------------- ---------
            Assert.Equal(haplotypeSegments.HaplotypeA.Count, 1);
            Assert.Equal(haplotypeSegments.HaplotypeB.Count, 3);
        }

        [Fact]
        public void TestSplitSeveralCommonCNVOverlapsCanvasCNV()
        {
            // scenario: Canvas segment spans more than one common segment
            // canvasSegment:   ------------------------------------------------
            // commonSegment:            ------------     -------------------
            var counts = new List<SampleGenomicBin>
            {
                new SampleGenomicBin("chr1", 100000, 100001, 100),
                new SampleGenomicBin("chr1", 150000, 150001, 90),
                new SampleGenomicBin("chr1", 200000, 200001, 110),
                new SampleGenomicBin("chr1", 250000, 250001, 100),
                new SampleGenomicBin("chr1", 300000, 300001, 95),
                new SampleGenomicBin("chr1", 350000, 350001, 105),
                new SampleGenomicBin("chr1", 400000, 400001, 105),
                new SampleGenomicBin("chr1", 450000, 450001, 105),
                new SampleGenomicBin("chr1", 500000, 500001, 105)
            };
            var canvasSegments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 100000, 500002, counts),
            };
            var commonSegments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 200000, 250001, counts.Skip(2).Take(2).ToList()),
                new CanvasSegment("chr1", 400000, 450001, counts.Skip(4).Take(2).ToList()),

            };
            var canvasSegmentsIndex = 0;
            var commonSegmentsIndex = 0;
            const int defaultReadCountsThreshold = 4;
            var haplotypeSegments = CanvasSegment.SplitCanvasSegments(canvasSegments, commonSegments, defaultReadCountsThreshold, ref canvasSegmentsIndex, ref commonSegmentsIndex);

            // transform into "haplotype" segments
            // canvasSegment:   ------------------------------------------------
            // commonSegment:   ---------  ------------     -------------------
            Assert.Equal(haplotypeSegments.HaplotypeA.Count, 1);
            Assert.Equal(haplotypeSegments.HaplotypeB.Count, 3);
        }

        [Fact]
        public void TestSplitCommonCNVPartOverlapsCanvasCNV()
        {
            // scenario: Canvas segment part overlaps common segment and comes first
            // canvasSegment:   --------------
            // commonSegment:            ------------     
            var counts = new List<SampleGenomicBin>
            {
                new SampleGenomicBin("chr1", 100000, 100001, 100),
                new SampleGenomicBin("chr1", 150000, 150001,  90),
                new SampleGenomicBin("chr1", 200000, 200001, 110),
                new SampleGenomicBin("chr1", 250000, 250001, 100),
                new SampleGenomicBin("chr1", 300000, 300001,  95),
                new SampleGenomicBin("chr1", 350000, 350001, 105),
                new SampleGenomicBin("chr1", 400000, 400001, 105),
                new SampleGenomicBin("chr1", 450000, 450001, 105),
                new SampleGenomicBin("chr1", 500000, 500001, 105)
            };
            var canvasSegments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 100000, 250001, counts.Take(4).ToList()),
                new CanvasSegment("chr1", 300000, 500001, counts.Skip(4).Take(5).ToList())
            };
            var commonSegments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 200000, 350001, counts.Skip(2).Take(4).ToList())
            };
            var canvasSegmentsIndex = 0;
            var commonSegmentsIndex = 0;
            const int defaultReadCountsThreshold = 4;
            var haplotypeSegments = CanvasSegment.SplitCanvasSegments(canvasSegments, commonSegments, defaultReadCountsThreshold, ref canvasSegmentsIndex, ref commonSegmentsIndex);

            // transform into "haplotype" segments
            // canvasSegment:   --------------  --------
            // commonSegment:   ---------   ------------ 
            Assert.Equal(haplotypeSegments.HaplotypeA.Count, 2);
            Assert.Equal(haplotypeSegments.HaplotypeB.Count, 2);
        }

        [Fact]
        public void TestSplitCommonCNVPartOverlapsCanvasCNVWithSameEndCoords()
        {
            // scenario: Canvas segment part overlaps common segment and comes first
            // canvasSegment:   ---------------------
            // commonSegment:            ------------     
            var counts = new List<SampleGenomicBin>
            {
                new SampleGenomicBin("chr1", 100000, 100001, 100),
                new SampleGenomicBin("chr1", 150000, 150001,  90),
                new SampleGenomicBin("chr1", 200000, 200001, 110),
                new SampleGenomicBin("chr1", 250000, 250001, 100),
                new SampleGenomicBin("chr1", 300000, 300001,  95),
                new SampleGenomicBin("chr1", 350000, 350001, 105),
                new SampleGenomicBin("chr1", 400000, 400001, 105),
                new SampleGenomicBin("chr1", 450000, 450001, 105),
                new SampleGenomicBin("chr1", 500000, 500001, 105)
            };
            var canvasSegments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 100000, 500001, counts),
            };
            var commonSegments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 300000, 500001, counts.Skip(4).Take(5).ToList())
            };
            var canvasSegmentsIndex = 0;
            var commonSegmentsIndex = 0;
            const int defaultReadCountsThreshold = 4;
            var haplotypeSegments = CanvasSegment.SplitCanvasSegments(canvasSegments, commonSegments, defaultReadCountsThreshold, ref canvasSegmentsIndex, ref commonSegmentsIndex);

            // transform into "haplotype" segments
            // canvasSegment:   ----------------------
            // commonSegment:   ---------   ---------- 
            Assert.Equal(haplotypeSegments.HaplotypeA.Count, 1);
            Assert.Equal(haplotypeSegments.HaplotypeB.Count, 2);
        }

        [Fact]
        public void TestSplitCommonCNVOverlapsSeveralCanvasCNV()
        {
            // scenario: common segment spans more than one Canvas segment
            // canvasSegment:   -------------------      --------
            // commonSegment:            ---------------------------------
            var counts = new List<SampleGenomicBin>
            {
                new SampleGenomicBin("chr1", 100000, 100001, 100),
                new SampleGenomicBin("chr1", 150000, 150001,  90),
                new SampleGenomicBin("chr1", 200000, 200001, 110),
                new SampleGenomicBin("chr1", 250000, 250001, 100),
                new SampleGenomicBin("chr1", 300000, 300001,  95),
                new SampleGenomicBin("chr1", 350000, 350001, 105),
                new SampleGenomicBin("chr1", 400000, 400001, 105),
                new SampleGenomicBin("chr1", 450000, 450001, 105),
                new SampleGenomicBin("chr1", 500000, 500001, 105)
            };
            var canvasSegments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 150000, 250001, counts.Skip(1).Take(2).ToList()),
                new CanvasSegment("chr1", 400000, 450001, counts.Skip(4).Take(2).ToList()),
            };
            var commonSegments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 200000, 500001, counts.Skip(2).Take(7).ToList())
            };
            var canvasSegmentsIndex = 0;
            var commonSegmentsIndex = 0;
            const int defaultReadCountsThreshold = 4;
            var haplotypeSegments = CanvasSegment.SplitCanvasSegments(canvasSegments, commonSegments, defaultReadCountsThreshold, ref canvasSegmentsIndex, ref commonSegmentsIndex);

            // transform into "haplotype" segments
            // canvasSegment:   -------------------      --------
            // commonSegment:   -------- ---------------------------------
            Assert.Equal(haplotypeSegments.HaplotypeA.Count, 2);
            Assert.Equal(haplotypeSegments.HaplotypeB.Count, 2);
        }
        
        [Fact]
        public void TestSplitCommonCNVPartOverlapsCanvasCNVEndComesFirst()
        {
            // scenario: Canvas segment part overlaps common segment and comes first
            // canvasSegment:             --------------
            // commonSegment:    ------------     
            var counts = new List<SampleGenomicBin>
            {
                new SampleGenomicBin("chr1", 100000, 100001, 100),
                new SampleGenomicBin("chr1", 150000, 150001,  90),
                new SampleGenomicBin("chr1", 200000, 200001, 110),
                new SampleGenomicBin("chr1", 250000, 250001, 100),
                new SampleGenomicBin("chr1", 300000, 300001,  95),
                new SampleGenomicBin("chr1", 350000, 350001, 105),
                new SampleGenomicBin("chr1", 400000, 400001, 105),
                new SampleGenomicBin("chr1", 450000, 450001, 105),
                new SampleGenomicBin("chr1", 500000, 500001, 105)
            };
            var canvasSegments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 300000, 450001, counts.Skip(4).Take(4).ToList())
            };
            var commonSegments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 200000, 350001, counts.Skip(2).Take(4).ToList())
            };
            var canvasSegmentsIndex = 0;
            var commonSegmentsIndex = 0;
            const int defaultReadCountsThreshold = 4;
            var haplotypeSegments = CanvasSegment.SplitCanvasSegments(canvasSegments, commonSegments, defaultReadCountsThreshold, ref canvasSegmentsIndex, ref commonSegmentsIndex);

            // transform into "haplotype" segments
            // canvasSegment:             --------------
            // commonSegment:    ------------  ---------
            Assert.Equal(haplotypeSegments.HaplotypeA.Count, 1);
            Assert.Equal(haplotypeSegments.HaplotypeB.Count, 2);
        }

        [Fact]
        public void TestMergeCommonCnvSegments()
        {
            const int alleleThreshold = 2;
            var counts = new List<SampleGenomicBin>
            {
                new SampleGenomicBin("chr1", 100000, 100001, 100),
                new SampleGenomicBin("chr1", 150000, 150001,  90),
                new SampleGenomicBin("chr1", 200000, 200001, 110),
                new SampleGenomicBin("chr1", 250000, 250001, 100),
                new SampleGenomicBin("chr1", 300000, 300001,  95),
                new SampleGenomicBin("chr1", 350000, 350001, 105),
                new SampleGenomicBin("chr1", 400000, 400001, 105),
                new SampleGenomicBin("chr1", 450000, 450001, 105),
                new SampleGenomicBin("chr1", 500000, 500001, 105)
            };

            // Canvas segment comes before common segment and does not overlap it
            var slice = counts.Skip(1).Take(3).ToList();
            var canvasSegments = new List<CanvasSegment> { new CanvasSegment(slice.First().GenomicBin.Chromosome, slice.First().Start, slice.Last().Stop, slice) };
            slice = counts.Skip(4).Take(2).ToList();
            var commonSegments = new List<CanvasSegment> { new CanvasSegment(slice.First().GenomicBin.Chromosome, slice.First().Start, slice.Last().Stop, slice) };
            var segmentHaplotypeses = CanvasSegment.MergeCommonCnvSegments(canvasSegments, commonSegments, alleleThreshold);
            Assert.Equal(2, segmentHaplotypeses.Count);
            Assert.Equal(canvasSegments, segmentHaplotypeses.First().HaplotypeA);
            Assert.Equal(null, segmentHaplotypeses.First().HaplotypeB);
            Assert.Equal(null, segmentHaplotypeses.Last().HaplotypeA);
            Assert.Equal(commonSegments, segmentHaplotypeses.Last().HaplotypeB);

            // common segment comes before Canvas segment and does not overlap it
            slice = counts.Skip(1).Take(3).ToList();
            commonSegments = new List<CanvasSegment> { new CanvasSegment(slice.First().GenomicBin.Chromosome, slice.First().Start, slice.Last().Stop, slice) };
            slice = counts.Skip(4).Take(2).ToList();
            canvasSegments = new List<CanvasSegment> { new CanvasSegment(slice.First().GenomicBin.Chromosome, slice.First().Start, slice.Last().Stop, slice) };
            segmentHaplotypeses = CanvasSegment.MergeCommonCnvSegments(canvasSegments, commonSegments, alleleThreshold);
            Assert.Equal(2, segmentHaplotypeses.Count);
            Assert.Equal(null, segmentHaplotypeses.First().HaplotypeA);
            Assert.Equal(commonSegments, segmentHaplotypeses.First().HaplotypeB);
            Assert.Equal(canvasSegments, segmentHaplotypeses.Last().HaplotypeA);
            Assert.Equal(null, segmentHaplotypeses.Last().HaplotypeB);

            // Canvas segment and common segment have the same coordinates
            slice = counts.Skip(1).Take(3).ToList();
            commonSegments = new List<CanvasSegment> { new CanvasSegment(slice.First().GenomicBin.Chromosome, slice.First().Start, slice.Last().Stop, slice) };
            slice = counts.Skip(1).Take(3).ToList();
            canvasSegments = new List<CanvasSegment> { new CanvasSegment(slice.First().GenomicBin.Chromosome, slice.First().Start, slice.Last().Stop, slice) };
            segmentHaplotypeses = CanvasSegment.MergeCommonCnvSegments(canvasSegments, commonSegments, alleleThreshold);
            Assert.Equal(1, segmentHaplotypeses.Count);
            Assert.Equal(null, segmentHaplotypeses.First().HaplotypeA);
            Assert.Equal(commonSegments, segmentHaplotypeses.First().HaplotypeB);

        }
    }
}