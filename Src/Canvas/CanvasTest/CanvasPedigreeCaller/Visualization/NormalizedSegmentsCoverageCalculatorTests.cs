using System.Collections.Generic;
using System.Linq;
using CanvasCommon;
using CanvasPedigreeCaller.Visualization;
using Isas.SequencingFiles.Bed;
using Xunit;

namespace CanvasTest
{
    public class NormalizedSegmentsCoverageCalculatorTests
    {
        [Fact]
        public void NormalizedSegmentsCoverageCalculator_NoBins_ReturnsNoBedGraphEntries()
        {
            var calculator = new NormalizedSegmentsCoverageCalculator();
            var segments = Enumerable.Empty<CanvasSegment>().ToList();

            var results = calculator.Calculate(segments);

            Assert.Empty(results);
        }

   
        [Fact]
        public void NormalizedSegmentsCoverageCalculator_OneSegmentOneBinCopyNumberZero_ReturnsBedGraphEntryWithZeroCoverage()
        {
            var calculator = new NormalizedSegmentsCoverageCalculator();
            var segments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 0, 1, new List<SampleGenomicBin>
                {
                    new SampleGenomicBin("chr1", 0, 1, 3f)
                })
                {
                    CopyNumber = 0,
                    Filter = CanvasFilter.PassFilter
                }
            };

            var result = calculator.Calculate(segments).Single();

            Assert.Equal("chr1", result.Chromosome);
            Assert.Equal(0, result.Interval.Start);
            Assert.Equal(1, result.Interval.End);
            Assert.Equal(0, result.Value);
        }


        [Fact]
        public void NormalizedSegmentsCoverageCalculator_TestMedianCoverages()
        {
            var calculator = new NormalizedSegmentsCoverageCalculator();
            var segments = new List<CanvasSegment>();

            // Median of 10 = 10
            var segment1 = GetSegment("chr1", 20, 30);
            AddBin(segment1, "chr1", 20, 30, 10);
            segments.Add(segment1);

            // Median of 20,30,50 = 30
            var segment2 = GetSegment("chr1", 40, 70);
            AddBin(segment2, "chr1", 40, 50, 20);
            AddBin(segment2, "chr1", 51, 60, 30);
            AddBin(segment2, "chr1", 61, 70, 50);

            segments.Add(segment2);

            // Median of 60,80 = 70
            var segment3 = GetSegment("chr2", 20, 50);
            AddBin(segment3, "chr2", 20, 30, 60);
            AddBin(segment3, "chr2", 40, 50, 80);
            segments.Add(segment3);

            var segment4 = GetSegment("chr3", 20, 50);
            AddBin(segment4, "chr3", 20, 30, 0);
            AddBin(segment4, "chr3", 40, 50, 0);
            segments.Add(segment4);

            var entries = calculator.Calculate(segments, 1).ToList();
            Assert.Equal(4, entries.Count());
            CheckBedEntry(entries[0], 20, 30, "chr1", 10);
            CheckBedEntry(entries[1], 40, 70, "chr1", 30);
            CheckBedEntry(entries[2], 20, 50, "chr2", 70);
            CheckBedEntry(entries[3], 20, 50, "chr3", 0);

            entries = calculator.Calculate(segments, 0.5).ToList();
            Assert.Equal(4, entries.Count());
            CheckBedEntry(entries[0], 20, 30, "chr1", 5);
            CheckBedEntry(entries[1], 40, 70, "chr1", 15);
            CheckBedEntry(entries[2], 20, 50, "chr2", 35);
            CheckBedEntry(entries[3], 20, 50, "chr3", 0);

        }

        private void CheckBedEntry(BedGraphEntry actual, int expectedStart, int expectedEnd, string expectedChr,
            decimal expectedValue)
        {
            Assert.Equal(expectedStart, actual.Interval.Start);
            Assert.Equal(expectedEnd, actual.Interval.End);
            Assert.Equal(expectedChr, actual.Chromosome);
            Assert.Equal(expectedValue, actual.Value);

        }

        private CanvasSegment GetSegment(string chr, int start, int end)
        {
            return new CanvasSegment(chr, start, end,
                new List<SampleGenomicBin>());
        }

        private static void AddBin(CanvasSegment segment, string chr, int start, int end, float count)
        {
            segment.GenomicBins.Add(new SampleGenomicBin(chr, start, end, count));
        }
    }
}