using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using CanvasCommon;
using Illumina.Common;
using Isas.Framework.DataTypes;
using Isas.SequencingFiles;
using Isas.SequencingFiles.Vcf;
using JetBrains.Annotations;
using Xunit;

namespace CanvasTest
{
    public class ReferencePloidyTests
    {
        [Fact]
        public void GetReferencePloidyIntervals_EmptyVcf_ReferencePloidyIs2()
        {
            var ploidyVcfIntervals = Enumerable.Empty<ReferencePloidyInterval>();
            var queryInterval = new ReferenceInterval("chrX", new Interval(1, 2));
            var expectedPloidyIntervals = new ReferencePloidyInterval(queryInterval, 2).Yield();

            AssertEqualPloidy(ploidyVcfIntervals, queryInterval, expectedPloidyIntervals);
        }

        [Fact]
        public void GetReferencePloidyIntervals_VcfPloidy1AndSameQueryInterval_PloidyIs1()
        {
            var queryInterval = new ReferenceInterval("chrX", new Interval(1, 2));
            var ploidyVcfIntervals = new ReferencePloidyInterval(queryInterval, 1).Yield().ToList();

            var expectedPloidyIntervals = ploidyVcfIntervals;
            AssertEqualPloidy(ploidyVcfIntervals, queryInterval, expectedPloidyIntervals);
        }

        [Theory]
        [InlineData(true)]
        [InlineData(false)]
        public void GetReferencePloidyIntervals_VcfPloidy1AndPartialOverlapQueryInterval_PloidyIs1And2(bool symbolicAltAllele)
        {
            var vcfInterval = new ReferenceInterval("chrX", new Interval(1, 1));
            var ploidyVcfIntervals = new ReferencePloidyInterval(vcfInterval, 1).Yield().ToList();
            var queryInterval = new ReferenceInterval("chrX", new Interval(1, 2));

            var expectedPloidyIntervals = new[]
            {
                CreatePloidyInterval("chrX", new Interval(1, 1), 1),
                CreatePloidyInterval("chrX", new Interval(2, 2), 2)
            };
            AssertEqualPloidy(ploidyVcfIntervals, queryInterval, expectedPloidyIntervals, symbolicAltAllele);
        }

        [Fact]
        public void GetReferencePloidyIntervals_VcfAdjacentPloidy1_SinglePloidy1()
        {
            var ploidyVcfIntervals = new[]
            {
                CreatePloidyInterval("chrX", new Interval(1, 1), 1),
                CreatePloidyInterval("chrX", new Interval(2, 2), 1)
            };
            var queryInterval = new ReferenceInterval("chrX", new Interval(1, 2));

            var expectedPloidyIntervals = new[]
            {
                CreatePloidyInterval("chrX", new Interval(1, 2), 1)
            };
            AssertEqualPloidy(ploidyVcfIntervals, queryInterval, expectedPloidyIntervals);
        }

        [Fact]
        public void GetReferencePloidyIntervals_VcfWithOverlappingPloidy_ThrowsArgumentException()
        {
            var ploidyVcfIntervals = new[]
            {
                CreatePloidyInterval("chrX", new Interval(1, 1), 2),
                CreatePloidyInterval("chrX", new Interval(1, 2), 2)
            };

            var exception = Record.Exception(() => LoadReferencePloidy(ploidyVcfIntervals));

            Assert.IsType<ArgumentException>(exception);
        }

        [Fact]
        public void GetReferencePloidyIntervals_VcfMultiplePloidyAndLargeQuery_ReturnMultiplePloidies()
        {
            var ploidyVcfIntervals = new[]
            {
                CreatePloidyInterval("chrX", new Interval(2, 2), 1),
                CreatePloidyInterval("chrX", new Interval(4, 4), 3)
            };
            var queryInterval = new ReferenceInterval("chrX", new Interval(1, 5));

            var expectedPloidyIntervals = new[]
            {
                CreatePloidyInterval("chrX", new Interval(1, 1), 2),
                CreatePloidyInterval("chrX", new Interval(2, 2), 1),
                CreatePloidyInterval("chrX", new Interval(3, 3), 2),
                CreatePloidyInterval("chrX", new Interval(4, 4), 3),
                CreatePloidyInterval("chrX", new Interval(5, 5), 2)
            };
            AssertEqualPloidy(ploidyVcfIntervals, queryInterval, expectedPloidyIntervals);
        }

        [Fact]
        public void GetReferencePloidy_QueryWithNoPloidy1Region_Returns2()
        {
            var ploidyVcfIntervals = Enumerable.Empty<ReferencePloidyInterval>();

            var queryInterval = new ReferenceInterval("chrX", new Interval(1, 2));
            AssertReferencePloidy(ploidyVcfIntervals, queryInterval, 2);
        }

        [Fact]
        public void GetReferencePloidy_QueryContainedWithinPloidy1Region_Returns1()
        {
            var ploidyVcfIntervals = CreatePloidyInterval("chrX", new Interval(1, 4), 1).Yield();
            var queryInterval = new ReferenceInterval("chrX", new Interval(2, 3));

            AssertReferencePloidy(ploidyVcfIntervals, queryInterval, 1);
        }

        [Fact]
        public void GetReferencePloidy_QueryOverlapsMultiplePloidyRegions_ThrowsArgumentException()
        {
            var ploidyVcfIntervals = new[]
            {
                CreatePloidyInterval("chrX", new Interval(2, 2), 1)
            };
            var queryInterval = new ReferenceInterval("chrX", new Interval(1, 2));
            var referencePloidy = LoadReferencePloidy(ploidyVcfIntervals);

            var exception = Record.Exception(() => referencePloidy.GetReferencePloidy(queryInterval));

            Assert.IsType<ArgumentException>(exception);
        }

        private ReferencePloidyInterval CreatePloidyInterval(string chromosome, Interval interval, int ploidy)
        {
            return new ReferencePloidyInterval(new ReferenceInterval(chromosome, interval), ploidy);
        }

        private static ReferencePloidy LoadReferencePloidy(IEnumerable<ReferencePloidyInterval> ploidyVcfIntervals,
            bool useSymbolicAltAlleleInPloidyVcf = true)
        {
            var sampleId = new SampleId("sampleId");
            var vcf = GetVcfAsString(sampleId, ploidyVcfIntervals, useSymbolicAltAlleleInPloidyVcf);

            using (var textReader = new StringReader(vcf))
            using (var reader = new GzipOrTextReader(textReader))
            {
                return ReferencePloidy.Load(reader, sampleId);
            }
        }

        [AssertionMethod]
        private void AssertEqualPloidy(
            IEnumerable<ReferencePloidyInterval> ploidyVcfIntervals, ReferenceInterval queryInterval,
            IEnumerable<ReferencePloidyInterval> expectedPloidyIntervals, bool useSymbolicAltAlleleInPloidyVcf = true)
        {
            var referencePloidy = LoadReferencePloidy(ploidyVcfIntervals, useSymbolicAltAlleleInPloidyVcf);
            var referencePloidyIntervals = referencePloidy.GetReferencePloidyIntervals(queryInterval).ToList();
            Assert.Equal(expectedPloidyIntervals, referencePloidyIntervals, new EqualityComparer());
        }

        [AssertionMethod]
        private void AssertReferencePloidy(
            IEnumerable<ReferencePloidyInterval> ploidyVcfIntervals, ReferenceInterval queryInterval, int expectedPloidy)
        {
            var referencePloidy = LoadReferencePloidy(ploidyVcfIntervals);
            var ploidy = referencePloidy.GetReferencePloidy(queryInterval);
            Assert.Equal(expectedPloidy, ploidy);
        }

        private static string GetVcfAsString(SampleId sampleId, IEnumerable<ReferencePloidyInterval> intervals, bool symbolicAllele)
        {
            var vcf =
                "##fileformat=VCFv4.1\n" +
                $"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\t{sampleId}\tS3\n";
            foreach (var interval in intervals)
            {
                var position = symbolicAllele
                    ? interval.Interval.Interval.OneBasedStart - 1
                    : interval.Interval.Interval.OneBasedStart;
                var altAllele = symbolicAllele ? "<CNV>" : ".";

                vcf +=
                    $"{interval.Interval.Chromosome}\t{position}\t.\tN\t{altAllele}\t.\tPASS\tEND={interval.Interval.Interval.OneBasedEnd}\tCN\t.\t{interval.ReferencePloidy}\t.\n";
            }
            return vcf;
        }

        private sealed class EqualityComparer : IEqualityComparer<ReferencePloidyInterval>
        {
            public bool Equals(ReferencePloidyInterval x, ReferencePloidyInterval y)
            {
                if (ReferenceEquals(x, y)) return true;
                if (ReferenceEquals(x, null)) return false;
                if (ReferenceEquals(y, null)) return false;
                if (x.GetType() != y.GetType()) return false;
                return x.Interval.Equals(y.Interval) && x.ReferencePloidy == y.ReferencePloidy;
            }

            public int GetHashCode(ReferencePloidyInterval obj)
            {
                unchecked
                {
                    return (obj.Interval.GetHashCode() * 397) ^ obj.ReferencePloidy;
                }
            }
        }
    }
}