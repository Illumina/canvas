using System.Collections.Generic;
using System.IO;
using System.Linq;
using CanvasCommon;
using CanvasPedigreeCaller;
using CanvasPedigreeCaller.Visualization;
using Isas.SequencingFiles;
using Xunit;

namespace CanvasTest
{
    public class NormalizedCoverageWriterTests
    {
        [Fact]
        public void NormalizedCoverageCalculator_NoBins_ReturnsNoBedGraphEntries()
        {
            var calculator = new NormalizedCoverageCalculator();
            var segments = new List<CanvasSegment>();

            var results = calculator.Calculate(segments);

            Assert.Empty(results);
        }

        [Fact]
        public void NormalizedCoverageCalculator_OneSegmentOneBin_ReturnsBedGraphEntryWithSegmentCopyNumber()
        {
            var calculator = new NormalizedCoverageCalculator();
            var segments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 1, 2, new List<SampleGenomicBin>
                {
                    new SampleGenomicBin("chr1", 0, 1, 3f),
                })
                {
                    CopyNumber = 2
                }
            };

            var results = calculator.Calculate(segments).ToList();
            var result = results.Single();

            Assert.Equal("chr1", result.Chromosome);
            Assert.Equal(0, result.Interval.Start);
            Assert.Equal(2, results.Single().Value);
        }

        [Fact]
        public void NormalizedCoverageCalculator_NoNormalization2_WritesBinCounts()
        {
            var normalizedCoverageCalculator = new NormalizedCoverageCalculator();
            var segments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 1, 2, new List<SampleGenomicBin>
                {
                    new SampleGenomicBin("chr1", 0, 1, 2f),
                    new SampleGenomicBin("chr1", 1, 2, 3f),
                })
            };

            var coverageFile = new StringWriter();
            coverageFile.NewLine = "\n";
            using (coverageFile)
            using (var writer = new BgzipOrStreamWriter(coverageFile))
            {
                //normalizedCoverageWriter.Write(segments, writer, 2);
            }

            Assert.Equal("chr1\t0\t1\t1\nchr1\t1\t2\t1.5\n", coverageFile.ToString());
        }
    }
}