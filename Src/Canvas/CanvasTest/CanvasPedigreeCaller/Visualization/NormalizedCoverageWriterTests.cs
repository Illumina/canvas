using System;
using System.Collections.Generic;
using System.Linq;
using CanvasCommon;
using CanvasPedigreeCaller.Visualization;
using Xunit;

namespace CanvasTest
{
    public class NormalizedBinsCoverageWriterTests
    {
        [Fact]
        public void NormalizedBinsCoverageCalculator_NoSegments_ReturnsNoBedGraphEntries()
        {
            var calculator = new NormalizedBinsCoverageCalculator();
            var segments = Enumerable.Empty<CanvasSegment>().ToList();

            var results = calculator.Calculate(segments);

            Assert.Empty(results);
        }

        [Fact]
        public void NormalizedBinsCoverageCalculator_SegmentWithNoBins_ReturnsNoBedGraphEntries()
        {
            var calculator = new NormalizedBinsCoverageCalculator();
            var segment = new CanvasSegment("chr1", 100, 120, new List<SampleGenomicBin>());
            var segments = new List<CanvasSegment>() {segment};

            // Returns no bins because it was given no bins (works if normalization precomputed)
            var results = calculator.Calculate(segments, 1);
            Assert.Empty(results);

            // Throws exception if it tries to compute normalization factor but has nothing to do it with
            // If this is a reasonable scenario, might want to make a nicer exception than the unhandled AggregateException
            // But I suspect this will never actually happen. For now it's a purely theoretical edge case.
            Assert.Throws<AggregateException>(() => calculator.Calculate(segments));
        }

        [Fact]
        public void
            NormalizedBinsCoverageCalculator_OneSegmentOneBinCopyNumberZero_ReturnsBedGraphEntryWithZeroCoverage()
        {
            var calculator = new NormalizedBinsCoverageCalculator();
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
        public void NormalizedCoverageCalculator_OneSegmentOneBin_ReturnsBedGraphEntryWithSegmentCopyNumber()
        {
            var calculator = new NormalizedBinsCoverageCalculator();
            var segments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 0, 1, new List<SampleGenomicBin>
                {
                    new SampleGenomicBin("chr1", 0, 1, 3f)
                })
                {
                    CopyNumber = 2,
                    Filter = CanvasFilter.PassFilter
                }
            };

            var result = calculator.Calculate(segments).Single();

            Assert.Equal(2, result.Value);
        }

        [Fact]
        public void
            NormalizedCoverageCalculator_OneSegmentPassOneSegmentFiltered_ReturnsBedGraphEntriesNormalizedByPassingSegmentCopyNumber()
        {
            var calculator = new NormalizedBinsCoverageCalculator();
            var segments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 0, 1, new List<SampleGenomicBin>
                {
                    new SampleGenomicBin("chr1", 0, 1, 3f)
                })
                {
                    CopyNumber = 1,
                    Filter = CanvasFilter.PassFilter
                },
                new CanvasSegment("chr1", 1, 2, new List<SampleGenomicBin>
                {
                    new SampleGenomicBin("chr1", 1, 2, 6f)
                })
                {
                    CopyNumber = 10,
                    Filter = CanvasFilter.Create(new[] {"Filtered"})
                }
            };

            var results = calculator.Calculate(segments).ToList();

            Assert.Equal(2, results.Count);
            Assert.Equal(1, results[0].Value);
            Assert.Equal(2, results[1].Value);
        }

        [Fact]
        public void
            NormalizedBinsCoverageCalculator_TwoSegmentsPassingEqualWeighting_ReturnsBedGraphEntriesNormalizedByMedianCopyNumber()
        {
            var calculator = new NormalizedBinsCoverageCalculator();
            var segments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 0, 1, new List<SampleGenomicBin>
                {
                    new SampleGenomicBin("chr1", 0, 1, 4f)
                })
                {
                    CopyNumber = 3,
                    Filter = CanvasFilter.PassFilter
                },
                new CanvasSegment("chr1", 1, 2, new List<SampleGenomicBin>
                {
                    new SampleGenomicBin("chr1", 1, 2, 8f)
                })
                {
                    CopyNumber = 2,
                    Filter = CanvasFilter.PassFilter
                }
            };

            var results = calculator.Calculate(segments).ToList();

            Assert.Equal(2, results.Count);

            // normalization factor is average of 3/4 and 2/8 = 0.5 
            Assert.Equal(2, results[0].Value);
            Assert.Equal(4, results[1].Value);
        }

        [Fact]
        public void NormalizedBinsCoverageCalculator_PrecomputedNormalizationFactor()
        {
            var calculator = new NormalizedBinsCoverageCalculator();
            var segments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 0, 1, new List<SampleGenomicBin>
                {
                    new SampleGenomicBin("chr1", 0, 1, 4f)
                })
                {
                    CopyNumber = 3,
                    Filter = CanvasFilter.PassFilter
                },
                new CanvasSegment("chr1", 1, 2, new List<SampleGenomicBin>
                {
                    new SampleGenomicBin("chr1", 1, 2, 8f)
                })
                {
                    CopyNumber = 2,
                    Filter = CanvasFilter.PassFilter
                }
            };

            var normalizationFactor = 0.5;
            var results = calculator.Calculate(segments, normalizationFactor).ToList();

            Assert.Equal(2, results.Count);
            Assert.Equal(2, results[0].Value);
            Assert.Equal(4, results[1].Value);

            normalizationFactor = 1;
            results = calculator.Calculate(segments, normalizationFactor).ToList();
            Assert.Equal(2, results.Count);
            Assert.Equal(4, results[0].Value);
            Assert.Equal(8, results[1].Value);

            normalizationFactor = 0.25;
            results = calculator.Calculate(segments, normalizationFactor).ToList();
            Assert.Equal(2, results.Count);
            Assert.Equal(1, results[0].Value);
            Assert.Equal(2, results[1].Value);

            results = calculator.Calculate(new List<CanvasSegment>(), normalizationFactor).ToList();
            Assert.Empty(results);

        }

    }
}