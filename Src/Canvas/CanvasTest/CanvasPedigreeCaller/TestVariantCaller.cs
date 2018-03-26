using System.Collections.Generic;
using System.Linq;
using CanvasCommon;
using CanvasPedigreeCaller;
using Isas.Framework.DataTypes;
using Isas.Framework.DataTypes.Maps;
using Xunit;

namespace CanvasTest.CanvasPedigreeCaller
{
    public class TestVariantCaller
    {
        [Fact]
        public void TestCommonCnvAssignment_DeNovoVariants()
        {
            var bins = new List<SampleGenomicBin>
            {
                new SampleGenomicBin("chr1", 1, 2, 100),
                new SampleGenomicBin("chr1", 1, 2, 100),
                new SampleGenomicBin("chr1", 1, 2, 100)
            };
            var balleles = new Balleles(new List<Ballele> {new Ballele(5501, 30, 30)});
            var segmentParent1 = new CanvasSegment("chr1", 1, 2, bins, balleles) {CopyNumber = 2};

            bins = new List<SampleGenomicBin>
            {
                new SampleGenomicBin("chr1", 1, 2, 100),
                new SampleGenomicBin("chr1", 1, 2, 100),
                new SampleGenomicBin("chr1", 1, 2, 100)
            };
            balleles = new Balleles(new List<Ballele> {new Ballele(5501, 30, 30)});
            var segmentParent2 = new CanvasSegment("chr1", 1, 2, bins, balleles) {CopyNumber = 2};

            bins = new List<SampleGenomicBin>
            {
                new SampleGenomicBin("chr1", 1, 2, 0),
                new SampleGenomicBin("chr1", 1, 2, 0),
                new SampleGenomicBin("chr1", 1, 2, 0)
            };
            balleles = new Balleles(new List<Ballele> {new Ballele(5501, 0, 0)});
            var segmentProband = new CanvasSegment("chr1", 1, 2, bins, balleles) {CopyNumber = 0};

            var pedigreeSegments = new SampleMap<CanvasSegment>
            {
                {new SampleId("parent1"), segmentParent1},
                {new SampleId("parent2"), segmentParent2},
                {new SampleId("proband"), segmentProband}
            };

            var sampleMetricsParent1 = SampleMetrics.GetSampleInfo(new List<CanvasSegment> {segmentParent1},
                ploidyBedPath: null,
                numberOfTrimmedBins: 2, id: new SampleId("parent1"));
            var sampleMetricsParent2 = SampleMetrics.GetSampleInfo(new List<CanvasSegment> {segmentParent2},
                ploidyBedPath: null,
                numberOfTrimmedBins: 2, id: new SampleId("parent2"));
            var sampleMetricsProband = SampleMetrics.GetSampleInfo(new List<CanvasSegment> {segmentProband},
                ploidyBedPath: null,
                numberOfTrimmedBins: 2, id: new SampleId("proband"));

            var sampleMetrics = new SampleMap<SampleMetrics>
            {
                {new SampleId("parent1"), sampleMetricsParent1},
                {new SampleId("parent2"), sampleMetricsParent2},
                {new SampleId("proband"), sampleMetricsProband}
            };

            bool isCommonCnv = global::CanvasPedigreeCaller.CanvasPedigreeCaller.IsSharedCnv(pedigreeSegments, sampleMetrics,
                new List<SampleId> {new SampleId("parent1"), new SampleId("parent2")},
                new SampleId("proband"), maximumCopyNumber: 5);

            Assert.False(isCommonCnv);

            var pedigreeGenotypes = new SampleMap<Genotype>
            {
                {new SampleId("parent1"), Genotype.Create(new PhasedGenotype(1, 1))},
                {new SampleId("parent2"), Genotype.Create(new PhasedGenotype(1, 1))},
                {new SampleId("proband"), Genotype.Create(new PhasedGenotype(0, 1))}
            };

            isCommonCnv = global::CanvasPedigreeCaller.CanvasPedigreeCaller.IsSharedCnv(pedigreeGenotypes, pedigreeSegments,
                sampleMetrics, new List <SampleId> {new SampleId("parent1"), new SampleId("parent2")},
                new SampleId("proband"), maximumCopyNumber: 5);

            Assert.False(isCommonCnv);

            pedigreeGenotypes = new SampleMap<Genotype>
            {
                {new SampleId("parent1"), Genotype.Create(new PhasedGenotype(2, 1))},
                {new SampleId("parent2"), Genotype.Create(new PhasedGenotype(1, 1))},
                {new SampleId("proband"), Genotype.Create(new PhasedGenotype(0, 1))}
            };

            isCommonCnv = global::CanvasPedigreeCaller.CanvasPedigreeCaller.IsSharedCnv(pedigreeGenotypes, pedigreeSegments,
                sampleMetrics, new List<SampleId> { new SampleId("parent1"), new SampleId("parent2") },
                new SampleId("proband"), maximumCopyNumber: 5);

            Assert.False(isCommonCnv);

        }

        [Fact]
        public void TestCommonCnvAssignment_InheritedVariants()
        {
            var bins = new List<SampleGenomicBin>
            {
                new SampleGenomicBin("chr1", 1, 2, 100),
                new SampleGenomicBin("chr1", 1, 2, 100),
                new SampleGenomicBin("chr1", 1, 2, 100)
            };
            var balleles = new Balleles(new List<Ballele> {new Ballele(5501, 30, 30)});
            var segmentParent1 = new CanvasSegment("chr1", 1, 2, bins, balleles) {CopyNumber = 2};

            bins = new List<SampleGenomicBin>
            {
                new SampleGenomicBin("chr1", 1, 2, 50),
                new SampleGenomicBin("chr1", 1, 2, 50),
                new SampleGenomicBin("chr1", 1, 2, 50)
            };
            balleles = new Balleles(new List<Ballele> {new Ballele(5501, 0, 30)});
            var segmentParent2 = new CanvasSegment("chr1", 1, 2, bins, balleles) {CopyNumber = 1};

            bins = new List<SampleGenomicBin>
            {
                new SampleGenomicBin("chr1", 1, 2, 50),
                new SampleGenomicBin("chr1", 1, 2, 50),
                new SampleGenomicBin("chr1", 1, 2, 50)
            };
            balleles = new Balleles(new List<Ballele> {new Ballele(5501, 0, 30)});
            var segmentProband = new CanvasSegment("chr1", 1, 2, bins, balleles) {CopyNumber = 1};

            var pedigreeSegments = new SampleMap<CanvasSegment>
            {
                {new SampleId("parent1"), segmentParent1},
                {new SampleId("parent2"), segmentParent2},
                {new SampleId("proband"), segmentProband}
            };

            var sampleMetricsParent1 = SampleMetrics.GetSampleInfo(new List<CanvasSegment> {segmentParent1},
                ploidyBedPath: null,
                numberOfTrimmedBins: 2, id: new SampleId("parent1"));
            var sampleMetricsParent2 = SampleMetrics.GetSampleInfo(new List<CanvasSegment> {segmentParent2},
                ploidyBedPath: null,
                numberOfTrimmedBins: 2, id: new SampleId("parent1"));
            var sampleMetricsProband = SampleMetrics.GetSampleInfo(new List<CanvasSegment> {segmentProband},
                ploidyBedPath: null,
                numberOfTrimmedBins: 2, id: new SampleId("parent1"));

            var sampleMetrics = new SampleMap<SampleMetrics>
            {
                {new SampleId("parent1"), sampleMetricsParent1},
                {new SampleId("parent2"), sampleMetricsParent2},
                {new SampleId("proband"), sampleMetricsProband}
            };

            bool isCommonCnv = global::CanvasPedigreeCaller.CanvasPedigreeCaller.IsSharedCnv(pedigreeSegments, sampleMetrics,
                new List<SampleId> {new SampleId("parent1"), new SampleId("parent2")},
                new SampleId("proband"), maximumCopyNumber: 5);

            Assert.True(isCommonCnv);

            var pedigreeGenotypes = new SampleMap<Genotype>
            {
                {new SampleId("parent1"), Genotype.Create(new PhasedGenotype(1, 1))},
                {new SampleId("parent2"), Genotype.Create(new PhasedGenotype(0, 1))},
                {new SampleId("proband"), Genotype.Create(new PhasedGenotype(0, 1))}
            };

            isCommonCnv = global::CanvasPedigreeCaller.CanvasPedigreeCaller.IsSharedCnv(pedigreeSegments, sampleMetrics,
                new List<SampleId> { new SampleId("parent1"), new SampleId("parent2") },
                new SampleId("proband"), maximumCopyNumber: 5);

            Assert.True(isCommonCnv);

            pedigreeGenotypes = new SampleMap<Genotype>
            {
                {new SampleId("parent1"), Genotype.Create(1)},
                {new SampleId("parent2"), Genotype.Create(1)},
                {new SampleId("proband"), Genotype.Create(0)}
            };

            global::CanvasPedigreeCaller.CanvasPedigreeCaller.IsSharedCnv(pedigreeSegments, sampleMetrics,
                new List<SampleId> { new SampleId("parent1"), new SampleId("parent2") },
                new SampleId("proband"), maximumCopyNumber: 5);

            Assert.True(isCommonCnv);
        }


        [Fact]
        public void TestGetGtLogLikelihoodScore()
        {
            var copyNumberModelFactory = new HaplotypeCopyNumberModelFactory();
            var copyNumberModel = copyNumberModelFactory.CreateModel(numCnStates: 5, maxCoverage: 200,
                meanCoverage: 100, diploidAlleleMeanCounts: 50.0);
            var simulatedCn = 3;
            var gtModelCounts = PedigreeInfo.GeneratePhasedGenotype(numCnStates: 5).Where(gt => gt.TotalCopyNumber == simulatedCn)
                .Select(gt => gt.PhasedGenotype).ToList();
            var gtObservedCounts = new Balleles(new List<Ballele>
            {
                new Ballele(1, 1, 73),
                new Ballele(100, 2, 74),
                new Ballele(200, 1, 76),
                new Ballele(300, 0, 74),
                new Ballele(400, 1, 75),
                new Ballele(500, 2, 74)
            });
            var gt0_3 = new PhasedGenotype(3, 0);
            int? selectedGtState = null;
            // variant caller only calls MCC, only upper triangle of CN genotypes is selected - i.e. CNa=3,CNb=0 from [CNa=3,CNb=0,CNa=0,CNb=3]
            double logLikelihoodScore =
                VariantCaller.GetGtLogLikelihoodScore(gtObservedCounts, gtModelCounts, ref selectedGtState,
                    copyNumberModel);
            Assert.Equal(gtModelCounts.IndexOf(gt0_3), selectedGtState);

            gtObservedCounts = new Balleles(new List<Ballele>
            {
                new Ballele(1, 23, 53),
                new Ballele(100, 22, 54),
                new Ballele(200, 25, 46),
                new Ballele(300, 24, 50),
                new Ballele(400, 26, 51),
                new Ballele(500, 24, 51)
            });
            var gt1_2 = new PhasedGenotype(2, 1);
            selectedGtState = null;
            // variant caller only calls MCC, only upper triangle of CN genotypes is selected - i.e. CNa=3,CNb=0 from [CNa=3,CNb=0,CNa=0,CNb=3]
            logLikelihoodScore =
                VariantCaller.GetGtLogLikelihoodScore(gtObservedCounts, gtModelCounts, ref selectedGtState,
                    copyNumberModel);
            Assert.Equal(gtModelCounts.IndexOf(gt1_2), selectedGtState);

        }
    }
}