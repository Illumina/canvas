using System.Collections.Generic;
using System.IO;
using System.Linq;
using CanvasCommon;
using Isas.SequencingFiles;
using Xunit;

namespace CanvasTest
{
    public class NormalizedCoverageWriterTests
    {
        [Fact]
        public void NormalizedCoverageWriter_NoBins_WritesEmptyFile()
        {
            var normalizedCoverageWriter = new NormalizedCoverageWriter();
            var segments = new List<CanvasSegment>();

            var coverageFile = new StringWriter();
            using (coverageFile)
            using (var writer = new BgzipOrStreamWriter(coverageFile))
            {
                normalizedCoverageWriter.Write(segments, writer, 1);
            }

            Assert.Empty(coverageFile.ToString());
        }

        [Fact]
        public void NormalizedCoverageWriter_NoNormalization_WritesBinCounts()
        {
            var normalizedCoverageWriter = new NormalizedCoverageWriter();
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
                normalizedCoverageWriter.Write(segments, writer, 1);
            }

            Assert.Equal("chr1\t0\t1\t2\nchr1\t1\t2\t3\n", coverageFile.ToString());
        }

        [Fact]
        public void NormalizedCoverageWriter_NoNormalization2_WritesBinCounts()
        {
            var normalizedCoverageWriter = new NormalizedCoverageWriter();
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
                normalizedCoverageWriter.Write(segments, writer, 2);
            }

            Assert.Equal("chr1\t0\t1\t1\nchr1\t1\t2\t1.5\n", coverageFile.ToString());
        }
    }

    public class NormalizedCoverageWriter
    {
        public void Write(IEnumerable<CanvasSegment> segments, BgzipOrStreamWriter coverageFile, double normalizationFactor)
        {
            foreach (var segment in segments)
            {
                foreach (var bin in segment.GenomicBins)
                {
                    coverageFile.WriteLine($"{bin.GenomicBin.Chromosome}\t{bin.GenomicBin.Interval.Start}\t{bin.GenomicBin.Interval.End}\t{bin.Count / normalizationFactor}");
                }
            }
        }
    }
}