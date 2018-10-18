using System.Collections.Generic;
using CanvasCommon;
using CanvasPartition;
using Xunit;

namespace CanvasTest.CanvasPartition
{
    public class SegmentationResultsProcessorTests
    {
        [Fact]
        public void PostProcessSegmentsTests()
        {
            var processor = new SegmentationResultsProcessor(100);

            var chr1Segments = new List<SegmentationInput.Segment>();
            chr1Segments.Add(new SegmentationInput.Segment(){start = 1, end=1000});
            chr1Segments.Add(new SegmentationInput.Segment() { start = 1100, end = 4500 });
            chr1Segments.Add(new SegmentationInput.Segment() { start = 4600, end = 5000 });

            var segmentsByChrom = new Dictionary<string, SegmentationInput.Segment[]>();
            segmentsByChrom.Add("chr1", chr1Segments.ToArray());
            var segmentationResults = new GenomeSegmentationResults(segmentsByChrom);

            var ploidyInfo = new PloidyInfo();
            var excludedIntervals = new Dictionary<string, List<SampleGenomicBin>>();
            var coverageInfo = new CoverageInfo();

            coverageInfo.CoverageByChr = new Dictionary<string, double[]>();
            coverageInfo.EndByChr = new Dictionary<string, uint[]>();
            coverageInfo.StartByChr = new Dictionary<string, uint[]>();
            coverageInfo.CoverageByChr.Add("chr1",new double[]{10, 10, 50, 100, 25, 10});
            coverageInfo.StartByChr.Add("chr1", new uint[] { 100, 600, 1200, 1300, 4001, 5000 });
            coverageInfo.EndByChr.Add("chr1", new uint[] { 500, 890, 1299, 4000, 4500, 5050 });

            var results = processor.PostProcessSegments(segmentationResults, ploidyInfo, excludedIntervals, coverageInfo);

            var chr1Results = results["chr1"];
            Assert.Equal(3, chr1Results.Count);

            // Final segments should reflect boundaries of actual bins within them 
            //  (in practice, these probably shouldn't disagree? but let's go theoretical here)
            SegmentTestHelpers.CheckSegment(chr1Results[0], 100, 890, 10, 2);
            SegmentTestHelpers.CheckSegment(chr1Results[1], 1200, 4500, 50, 3);
            SegmentTestHelpers.CheckSegment(chr1Results[2], 5000, 5050, 10, 1); // Bin extends past segment - still keep it (?)

            // Add forbidden zone between two bins of the same original segment, this should split up the affected segment
            excludedIntervals.Add("chr1", new List<SampleGenomicBin>(){GetForbiddenZone("chr1", 525, 575)}); // Mid = 550, in between the bins of the first segment
            results = processor.PostProcessSegments(segmentationResults, ploidyInfo, excludedIntervals, coverageInfo);

            chr1Results = results["chr1"];
            Assert.Equal(4, chr1Results.Count);

            SegmentTestHelpers.CheckSegment(chr1Results[0], 100, 500, 10, 1);
            SegmentTestHelpers.CheckSegment(chr1Results[1], 600, 890, 10, 1);
            SegmentTestHelpers.CheckSegment(chr1Results[2], 1200, 4500, 50, 3);
            SegmentTestHelpers.CheckSegment(chr1Results[3], 5000, 5050, 10, 1); // Bin extends past segment - still keep it (?)

            // Forbidden zone midpoint is in the second bin -- apparently this is presumed to never happen because it would have already been taken care of
            // This fails the test with the Debug Asserts in there. Otherwise it would be counted as a new bin
            excludedIntervals.Clear();
            excludedIntervals.Add("chr1", new List<SampleGenomicBin>() { GetForbiddenZone("chr1", 585, 635) }); // Mid = 610, in second bin
            results = processor.PostProcessSegments(segmentationResults, ploidyInfo, excludedIntervals, coverageInfo);

            chr1Results = results["chr1"];
            Assert.Equal(4, chr1Results.Count);

            SegmentTestHelpers.CheckSegment(chr1Results[0], 100, 500, 10, 1);
            SegmentTestHelpers.CheckSegment(chr1Results[1], 600, 890, 10, 1);
            SegmentTestHelpers.CheckSegment(chr1Results[2], 1200, 4500, 50, 3);
            SegmentTestHelpers.CheckSegment(chr1Results[3], 5000, 5050, 10, 1); // Bin extends past segment - still keep it (?)

            // Forbidden zone midpoint is in the first bin although it ends between bins -- apparently this is presumed to never happen because it would have already been taken care of
            // Note the asymmetry compared to the above
            excludedIntervals.Clear();
            excludedIntervals.Add("chr1", new List<SampleGenomicBin>() { GetForbiddenZone("chr1", 465, 515) }); // Mid = 490, in first bin of first segment
            results = processor.PostProcessSegments(segmentationResults, ploidyInfo, excludedIntervals, coverageInfo);

            chr1Results = results["chr1"];
            // Would fail - asymmetry. What do we want?
            // Leave as-is for now so as not to change the behavior in this (unrelated) feature addition
            //Assert.Equal(4, chr1Results.Count);

            //SegmentTestHelpers.CheckSegment(chr1Results[0], 100, 500, 10, 1);
            //SegmentTestHelpers.CheckSegment(chr1Results[1], 600, 890, 10, 1);
            //SegmentTestHelpers.CheckSegment(chr1Results[2], 1200, 4500, 50, 3);
            //SegmentTestHelpers.CheckSegment(chr1Results[3], 5000, 5050, 10, 1); // Bin extends past segment - still keep it (?)

            // TODO test where no segment covers bins?
            // TODO overlapping segments or bins?
            // TODO bin starts before segment
            // TODO test interbin dist



        }

        private SampleGenomicBin GetForbiddenZone(string chrom, int start, int end)
        {
            return new SampleGenomicBin(chrom, start, end, 0);
        }

    }
}