using System.Collections.Generic;
using System.Linq;
using CanvasCommon;
using CanvasPedigreeCaller.Visualization;
using Xunit;

namespace CanvasTest
{
    public class NormalizedCoverageWriterTests
    {
        [Fact]
        public void NormalizedCoverageCalculator_NoBins_ReturnsNoBedGraphEntries()
        {
            var calculator = new NormalizedCoverageCalculator();
            var segments = Enumerable.Empty<CanvasSegment>().ToList();

            var results = calculator.Calculate(segments);

            Assert.Empty(results);
        }

        [Fact]
        public void NormalizedCoverageCalculator_OneSegmentOneBinCopyNumberZero_ReturnsBedGraphEntryWithZeroCoverage()
        {
            var calculator = new NormalizedCoverageCalculator();
            var segments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 0, 1, new List<SampleGenomicBin>
                {
                    new SampleGenomicBin("chr1", 0, 1, 3f)
                })
                {
                    CopyNumber = 0
                }
            };

            var result = calculator.Calculate(segments).Single();

            Assert.Equal("chr1", result.Chromosome);
            Assert.Equal(0, result.Interval.Start);
            Assert.Equal(1, result.Interval.End);
            Assert.Equal(0, result.Value);
        }

        [Fact]
        public void NormalizedCoverageCalculator_OneSegmentOneBin_ReturnsBedGraphEntryWithSegmentCopyNumber()
        {
            var calculator = new NormalizedCoverageCalculator();
            var segments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 0, 1, new List<SampleGenomicBin>
                {
                    new SampleGenomicBin("chr1", 0, 1, 3f)
                })
                {
                    CopyNumber = 2
                }
            };

            var result = calculator.Calculate(segments).Single();

            Assert.Equal(2, result.Value);
        }

        [Fact]
        public void NormalizedCoverageCalculator_OneSegmentPassOneSegmentFiltered_ReturnsBedGraphEntriesNormalizedByPassingSegmentCopyNumber()
        {
            var calculator = new NormalizedCoverageCalculator();
            var segments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 0, 1, new List<SampleGenomicBin>
                {
                    new SampleGenomicBin("chr1",0, 1, 3f)
                })
                {
                    CopyNumber = 1,
                    Filter = "PASS"
                },
                new CanvasSegment("chr1", 1, 2, new List<SampleGenomicBin>
                {
                    new SampleGenomicBin("chr1", 1, 2, 6f)
                })
                {
                    CopyNumber = 10,
                    Filter = "Filtered"
                }
            };

            var results = calculator.Calculate(segments).ToList();

            Assert.Equal(2, results.Count);
            Assert.Equal(1, results[0].Value);
            Assert.Equal(2, results[1].Value);
        }

        [Fact]
        public void NormalizedCoverageCalculator_TwoSegmentsPassingEqualWeighting_ReturnsBedGraphEntriesNormalizedByMedianCopyNumber()
        {
            var calculator = new NormalizedCoverageCalculator();
            var segments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 0, 1, new List<SampleGenomicBin>
                {
                    new SampleGenomicBin("chr1", 0, 1, 4f)
                })
                {
                    CopyNumber = 3,
                    Filter = "PASS"
                },
                new CanvasSegment("chr1", 1, 2, new List<SampleGenomicBin>
                {
                    new SampleGenomicBin("chr1", 1, 2, 8f)
                })
                {
                    CopyNumber = 2,
                    Filter = "PASS"
                }
            };

            var results = calculator.Calculate(segments).ToList();

            Assert.Equal(2, results.Count);

            // normalization factor is average of 3/4 and 2/8 = 0.5 
            Assert.Equal(2, results[0].Value);
            Assert.Equal(4, results[1].Value);
        }
    }
}