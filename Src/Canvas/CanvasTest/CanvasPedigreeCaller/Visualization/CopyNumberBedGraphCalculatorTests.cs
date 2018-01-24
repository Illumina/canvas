using System.Collections.Generic;
using System.Linq;
using CanvasCommon;
using CanvasPedigreeCaller.Visualization;
using Illumina.Common;
using Isas.SequencingFiles;
using Xunit;

namespace CanvasTest
{
    public class CopyNumberBedGraphCalculatorTests
    {
        [Fact]
        public void NoSegments_ReturnsNoBedGraphEntries()
        {
            var calculator = new CopyNumberBedGraphCalculator();
            var segments = Enumerable.Empty<CanvasSegment>().ToList();
            var ploidyInfo = new PloidyInfo();

            var results = calculator.Calculate(segments, ploidyInfo);

            Assert.Empty(results);
        }

        [Fact]
        public void FiltersNonPassSegments()
        {
            var calculator = new CopyNumberBedGraphCalculator();
            var segments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 0, 1, new List<SampleGenomicBin>
                {
                    new SampleGenomicBin("chr1", 0, 1, 3f)
                })
                {
                    CopyNumber = 0,
                    Filter = CanvasFilter.Create("NonPass".Yield())
                }
            };
            var ploidyInfo = new PloidyInfo();

            var results = calculator.Calculate(segments, ploidyInfo);

            Assert.Empty(results);
        }

        [Fact]
        public void VariantCopyNumber_ReturnsCopyNumber()
        {
            var calculator = new CopyNumberBedGraphCalculator();
            var segments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 0, 1, new List<SampleGenomicBin>
                {
                    new SampleGenomicBin("chr1", 0, 1, 3f)
                })
                {
                    CopyNumber = 1,
                    Filter = CanvasFilter.PassFilter
                }
            };
            var ploidyInfo = new PloidyInfo();

            var results = calculator.Calculate(segments, ploidyInfo).ToList();

            Assert.Equal("chr1", results.First().Chromosome);
            Assert.Equal(new BedInterval(0, 1), results.First().Interval);
            Assert.Equal(1m, results.First().Value);
        }

        [Fact]
        public void ReferenceCopyNumber_IsExcluded()
        {
            var calculator = new CopyNumberBedGraphCalculator();
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
            var ploidyInfo = new PloidyInfo();

            var results = calculator.Calculate(segments, ploidyInfo).ToList();

            Assert.Empty(results);
        }

        [Fact]
        public void ReferenceCopyNumberByPloidy_IsExcluded()
        {
            var calculator = new CopyNumberBedGraphCalculator();
            var segments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 0, 1, new List<SampleGenomicBin>
                {
                    new SampleGenomicBin("chr1", 0, 1, 3f)
                })
                {
                    CopyNumber = 1,
                    Filter = CanvasFilter.PassFilter
                }
            };
            var ploidyInfo = new PloidyInfo { PloidyByChromosome = ("chr1", new PloidyInterval("chr1") { Start = 0, End = 1, Ploidy = 1 }.Yield().ToList()).Yield().ToDictionary() };
            
            var results = calculator.Calculate(segments, ploidyInfo).ToList();

            Assert.Empty(results);
        }

        [Fact]
        public void LOH_IsIncluded()
        {
            var calculator = new CopyNumberBedGraphCalculator();
            var segments = new List<CanvasSegment>
            {
                new CanvasSegment("chr1", 0, 1, new List<SampleGenomicBin>
                {
                    new SampleGenomicBin("chr1", 0, 1, 3f)
                })
                {
                    CopyNumber = 2,
                    MajorChromosomeCount = 2,
                    Filter = CanvasFilter.PassFilter
                }
            };
            var ploidyInfo = new PloidyInfo();

            var results = calculator.Calculate(segments, ploidyInfo).ToList();
            
            Assert.Equal(2m, results.First().Value);
        }
    }
}