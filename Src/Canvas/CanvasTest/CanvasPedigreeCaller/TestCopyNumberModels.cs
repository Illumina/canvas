using System.Collections.Generic;
using CanvasPedigreeCaller;
using CanvasCommon;
using Xunit;

namespace CanvasTest.CanvasPedigreeCaller
{
    public class TestCopyNumberModels
    {
        [Fact]
        public void copyNumberModelTester_TotalCopyNumberOnly()
        {
            var copyNumberModelFactory = new CopyNumberModelFactory();
            var copyNumberModel = copyNumberModelFactory.CreateModel(numCnStates: 5, maxCoverage: 200, meanCoverage: 100, diploidAlleleMeanCounts: 50.0);
            double diploidLikelihood = copyNumberModel.GetTotalCopyNumberLikelihoods(100.0, Genotype.Create(2));
            double haploidLikelihood = copyNumberModel.GetTotalCopyNumberLikelihoods(100.0, Genotype.Create(1));
            Assert.True(diploidLikelihood > haploidLikelihood);
        }

        [Fact]
        public void haplotypeCopyNumberModelTester_PhasedGenotype_LossOfHeterozygosity()
        {
            var copyNumberModelFactory = new HaplotypeCopyNumberModelFactory();
            var copyNumberModel = copyNumberModelFactory.CreateModel(numCnStates: 5, maxCoverage: 200, meanCoverage: 100, diploidAlleleMeanCounts: 50.0);
            var bAlleles = new Balleles(new List<Ballele>
            {
                new Ballele(1, 50, 1),
                new Ballele(100, 25, 24),
                new Ballele(200, 23, 27),
                new Ballele(300, 25, 24),
                new Ballele(400, 1, 50),
                new Ballele(500, 25, 25)
            });
            double diploidHet = copyNumberModel.GetGenotypeLogLikelihood(bAlleles, new PhasedGenotype(1, 1));
            double lohB = copyNumberModel.GetGenotypeLogLikelihood(bAlleles, new PhasedGenotype(0, 2));
            double lohA = copyNumberModel.GetGenotypeLogLikelihood(bAlleles, new PhasedGenotype(2, 0));
            Assert.True(diploidHet > lohB);
            Assert.True(diploidHet > lohA);

            var bAllelesLohWithNoise = new Balleles(new List<Ballele>
            {
                new Ballele(1, 53, 1),
                new Ballele(100, 50, 1),
                new Ballele(200, 47, 2),
                new Ballele(300, 46, 0),
                new Ballele(400, 48, 2),
                new Ballele(500, 53, 0)
            });
            diploidHet = copyNumberModel.GetGenotypeLogLikelihood(bAllelesLohWithNoise, new PhasedGenotype(1, 1));
            lohB = copyNumberModel.GetGenotypeLogLikelihood(bAllelesLohWithNoise, new PhasedGenotype(0, 2));
            lohA = copyNumberModel.GetGenotypeLogLikelihood(bAllelesLohWithNoise, new PhasedGenotype(2, 0));
            Assert.True(diploidHet > lohB);
            Assert.True(diploidHet < lohA);
        }
    }
}