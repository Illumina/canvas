using System.Collections.Generic;
using System.Linq;
using Illumina.Common;
using JetBrains.Annotations;
using Xunit;
using GenomeSegmentationResults = CanvasPartition.GenomeSegmentationResults;
using Segment = CanvasPartition.SegmentationInput.Segment;

namespace CanvasTest.CanvasPartition
{
    public class GenomeSegmentationResultsTests
    {
        [Fact]
        public void SingleSamplesReturnsSameSegments()
        {
            var sampleSegmentationResults = new GenomeSegmentationResults(new Dictionary<string, Segment[]>
            {
                ["chr1"] = new[] { new Segment { start = 1, end = 200 } }
            });

            var segmentationResults =
                GenomeSegmentationResults.SplitOverlappingSegments(sampleSegmentationResults.Yield().ToList());

            AssertEqualSegmentations(sampleSegmentationResults, segmentationResults);
        }

        [Fact]
        public void IdenticalSegmentsReturnsOneCopy()
        {
            var sample1SegmentationResults = new GenomeSegmentationResults(new Dictionary<string, Segment[]>
            {
                ["chr1"] = new[] { new Segment { start = 1, end = 200 } }
            });
            var sample2SegmentationResults = new GenomeSegmentationResults(new Dictionary<string, Segment[]>
            {
                ["chr1"] = new[] { new Segment { start = 1, end = 200 } }
            });

            var segmentationResults =
                GenomeSegmentationResults.SplitOverlappingSegments(sample1SegmentationResults.Yield().Concat(sample2SegmentationResults).ToList());

            AssertEqualSegmentations(sample1SegmentationResults, sample2SegmentationResults);
            AssertEqualSegmentations(sample1SegmentationResults, segmentationResults);
        }

        [Fact]
        public void RecurringBoundariesAreUniquified()
        {
            var sample1SegmentationResults = new GenomeSegmentationResults(new Dictionary<string, Segment[]>
            {
                ["chr1"] = new[] { new Segment { start = 1, end = 300 } }
            });
            var sample2SegmentationResults = new GenomeSegmentationResults(new Dictionary<string, Segment[]>
            {
                ["chr1"] = new[] { new Segment { start = 1, end = 200 }, new Segment { start = 200, end = 300 } }
            });
            var sample3SegmentationResults = new GenomeSegmentationResults(new Dictionary<string, Segment[]>
            {
                ["chr1"] = new[] { new Segment { start = 1, end = 200 }, new Segment { start = 200, end = 300 } }
            });

            var segmentationResults =
                GenomeSegmentationResults.SplitOverlappingSegments(sample1SegmentationResults.Yield().Concat(sample2SegmentationResults).ToList());

            var expectedSegmentationResults = new GenomeSegmentationResults(new Dictionary<string, Segment[]>
            {
                ["chr1"] = new[] { new Segment { start = 1, end = 200 }, new Segment { start = 200, end = 300 } }
            });
            AssertEqualSegmentations(expectedSegmentationResults, segmentationResults);
        }

        [Fact]
        public void TwoPartialOverlappingSegments_ReturnsThreeSegments()
        {
            var sample1SegmentationResults = new GenomeSegmentationResults(new Dictionary<string, Segment[]>
            {
                ["chr1"] = new[] { new Segment { start = 0, end = 200 } }
            });
            var sample2SegmentationResults = new GenomeSegmentationResults(new Dictionary<string, Segment[]>
            {
                ["chr1"] = new[] { new Segment { start = 100, end = 300 } }
            });

            var segmentationResults =
                GenomeSegmentationResults.SplitOverlappingSegments(sample1SegmentationResults.Yield().Concat(sample2SegmentationResults).ToList());

            var expectedSegmentationResults = new GenomeSegmentationResults(new Dictionary<string, Segment[]>
            {
                ["chr1"] = new[] { new Segment { start = 0, end = 100 }, new Segment { start = 100, end = 200 }, new Segment { start = 200, end = 300 } }
            });
            AssertEqualSegmentations(expectedSegmentationResults, segmentationResults);
        }

        [Fact]
        public void TwoSegmentsSameStartPosition_ReturnsTwoSegments()
        {
            var sample1SegmentationResults = new GenomeSegmentationResults(new Dictionary<string, Segment[]>
            {
                ["chr1"] = new[] { new Segment { start = 0, end = 100 } }
            });
            var sample2SegmentationResults = new GenomeSegmentationResults(new Dictionary<string, Segment[]>
            {
                ["chr1"] = new[] { new Segment { start = 0, end = 200 } }
            });

            var segmentationResults =
                GenomeSegmentationResults.SplitOverlappingSegments(sample1SegmentationResults.Yield().Concat(sample2SegmentationResults).ToList());

            var expectedSegmentationResults = new GenomeSegmentationResults(new Dictionary<string, Segment[]>
            {
                ["chr1"] = new[] { new Segment { start = 0, end = 100 }, new Segment { start = 100, end = 200 } }
            });
            AssertEqualSegmentations(expectedSegmentationResults, segmentationResults);
        }

        [Fact]
        public void OneSegmentContainsAnotherSegment_ReturnsThreeSegments()
        {
            var sample1SegmentationResults = new GenomeSegmentationResults(new Dictionary<string, Segment[]>
            {
                ["chr1"] = new[] { new Segment { start = 0, end = 300 } }
            });
            var sample2SegmentationResults = new GenomeSegmentationResults(new Dictionary<string, Segment[]>
            {
                ["chr1"] = new[] { new Segment { start = 100, end = 200 } }
            });

            var segmentationResults =
                GenomeSegmentationResults.SplitOverlappingSegments(sample1SegmentationResults.Yield().Concat(sample2SegmentationResults).ToList());

            var expectedSegmentationResults = new GenomeSegmentationResults(new Dictionary<string, Segment[]>
            {
                ["chr1"] = new[] { new Segment { start = 0, end = 100 }, new Segment { start = 100, end = 200 }, new Segment { start = 200, end = 300 } }
            });
            AssertEqualSegmentations(expectedSegmentationResults, segmentationResults);
        }

        [Fact]
        public void OneSegmentOverlappingManySegments_ReturnsUnionOfBreakpoints()
        {
            var sample1SegmentationResults = new GenomeSegmentationResults(new Dictionary<string, Segment[]>
            {
                ["chr1"] = new[]
                {
                    new Segment { start = 0, end = 600 }
                }
            });
            var sample2SegmentationResults = new GenomeSegmentationResults(new Dictionary<string, Segment[]>
            {
                ["chr1"] = new[] 
                {
                    new Segment { start = 100, end = 200 },
                    new Segment { start = 200, end = 300 },
                    new Segment { start = 500, end = 700 },
                    new Segment { start = 800, end = 900 }
                }
            });

            var segmentationResults =
                GenomeSegmentationResults.SplitOverlappingSegments(sample1SegmentationResults.Yield().Concat(sample2SegmentationResults).ToList());

            var expectedSegmentationResults = new GenomeSegmentationResults(new Dictionary<string, Segment[]>
            {
                ["chr1"] = new[]
                {
                    new Segment { start = 0, end = 100 },
                    new Segment { start = 100, end = 200 },
                    new Segment { start = 200, end = 300 },
                    new Segment { start = 300, end = 500 },
                    new Segment { start = 500, end = 600 },
                    new Segment { start = 600, end = 700 },
                    new Segment { start = 800, end = 900 }
                }
            });
            AssertEqualSegmentations(expectedSegmentationResults, segmentationResults);
        }

        [Fact]
        public void OneSegmentEndsAtAnotherSegment_ReturnsTwoSegments()
        {
            var sample1SegmentationResults = new GenomeSegmentationResults(new Dictionary<string, Segment[]>
            {
                ["chr1"] = new[] { new Segment { start = 0, end = 100 } }
            });
            var sample2SegmentationResults = new GenomeSegmentationResults(new Dictionary<string, Segment[]>
            {
                ["chr1"] = new[] { new Segment { start = 100, end = 200 } }
            });

            var segmentationResults =
                GenomeSegmentationResults.SplitOverlappingSegments(sample1SegmentationResults.Yield().Concat(sample2SegmentationResults).ToList());

            var expectedSegmentationResults = new GenomeSegmentationResults(new Dictionary<string, Segment[]>
            {
                ["chr1"] = new[] { new Segment { start = 0, end = 100 }, new Segment { start = 100, end = 200 } }
            });
            AssertEqualSegmentations(expectedSegmentationResults, segmentationResults);
        }

        [Fact]
        public void CanHandleMultipleChromosomes()
        {
            var sample1SegmentationResults = new GenomeSegmentationResults(new Dictionary<string, Segment[]>
            {
                ["chr1"] = new[] { new Segment { start = 0, end = 100 } },
                ["chr2"] = new[] { new Segment { start = 300, end = 400 } }
            });
            var sample2SegmentationResults = new GenomeSegmentationResults(new Dictionary<string, Segment[]>
            {
                ["chr1"] = new[] { new Segment { start = 100, end = 200 } },
                ["chr2"] = new[] { new Segment { start = 500, end = 600 } }
            });

            var segmentationResults =
                GenomeSegmentationResults.SplitOverlappingSegments(sample1SegmentationResults.Yield().Concat(sample2SegmentationResults).ToList());

            var expectedSegmentationResults = new GenomeSegmentationResults(new Dictionary<string, Segment[]>
            {
                ["chr1"] = new[] { new Segment { start = 0, end = 100 }, new Segment { start = 100, end = 200 } },
                ["chr2"] = new[] { new Segment { start = 300, end = 400 }, new Segment { start = 500, end = 600 } }
            });
            AssertEqualSegmentations(expectedSegmentationResults, segmentationResults);
        }

        [AssertionMethod]
        private void AssertEqualSegmentations(GenomeSegmentationResults expected, GenomeSegmentationResults actual)
        {
            Assert.Equal(expected.SegmentByChr.Keys.ToHashSet(), actual.SegmentByChr.Keys.ToHashSet());
            foreach (var chr in expected.SegmentByChr.Keys)
            {
                IEqualityComparer<Segment> segmentComparer = new SegmentEqualityComparer();
                Assert.Equal(expected.SegmentByChr[chr], actual.SegmentByChr[chr], segmentComparer);
            }
        }

        private class SegmentEqualityComparer : IEqualityComparer<Segment>
        {
            public bool Equals(Segment x, Segment y)
            {
                return x.start == y.start && x.end == y.end;
            }

            public int GetHashCode(Segment obj)
            {
                return (obj.start, obj.end).GetHashCode();
            }
        }
    }
}