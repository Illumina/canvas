using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using CanvasCommon;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Isas.SequencingFiles;
using Isas.SequencingFiles.Bed;
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
                //normalizedCoverageWriter.Write(segments, writer, 1);
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
                //normalizedCoverageWriter.Write(segments, writer, 1);
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
                //normalizedCoverageWriter.Write(segments, writer, 2);
            }

            Assert.Equal("chr1\t0\t1\t1\nchr1\t1\t2\t1.5\n", coverageFile.ToString());
        }
    }

    public class NormalizedCoverageWriter
    {
        public void Write(IEnumerable<CanvasSegment> segments, IFileLocation coverageFile)
        {
            var normalizationFactor = ComputeWeightedNormalizationFactor(segments);
        }

        /// <summary>
        /// iterate through segments keeping a weighted sum of the normalization factor
        /// </summary>
        /// <param name="segments"></param>
        /// <returns></returns>
        private double ComputeWeightedNormalizationFactor(IEnumerable<CanvasSegment> segments)
        {
            var weightedSum = 0d;
            foreach (var segment in segments)
            {
                var normalizationFactor = ComputeNormalizationFactor(segment);
            }
            return weightedSum;
        }

        private double ComputeNormalizationFactor(CanvasSegment segment)
        {
            throw new NotSupportedException();
        }
    }

    public class BedGraphWriter
    {
        private readonly BgzipOrStreamWriter _writer;

        public BedGraphWriter(BgzipOrStreamWriter writer)
        {
            _writer = writer;
        }

        public void WriteLine(BedGraphEntry entry)
        {
            var bedGraphLine = GetBedGraphLine(entry);
            _writer.WriteLine(bedGraphLine);
        }

        private static string GetBedGraphLine(BedGraphEntry entry)
        {
            return $"{entry.Chromosome}\t{entry.Interval.Start}\t{entry.Interval.End}\t{entry.Value}";
        }
    }

    public static class BedGraphWriterExtensions
    {
        public static void WriteLines(this BedGraphWriter writer, IEnumerable<BedGraphEntry> entries)
        {
            entries.ForEach(writer.WriteLine);
        }

        public static void WriteLines(this IFileLocation location, IEnumerable<BedGraphEntry> entries)
        {
            using (var streamWriter = new BgzipOrStreamWriter(location.FullName))
            {
                var writer = new BedGraphWriter(streamWriter);
                writer.WriteLines(entries);
            }
        }
    }
}