using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CanvasCommon;
using CanvasPedigreeCaller.Visualization;
using Xunit;

namespace CanvasTest.CanvasPedigreeCaller
{
    public class TestGenotypes
    {
        [Fact]
        public void ContainsSharedAllelesTester_TotalCopyNumberOnly()
        {
            var genotype1 = Genotype.Create(3);
            var genotype2 = Genotype.Create(2);
            Assert.NotEqual(genotype1, genotype2);

            genotype2 = Genotype.Create(new PhasedGenotype(2,1));
            Assert.Equal(genotype1, genotype2);

            genotype2 = Genotype.Create(3);
            Assert.Equal(genotype1, genotype2);
        }

        [Fact]
        public void ContainsSharedAllelesTester_PhasedGenotype()
        {
            var genotype1 = Genotype.Create(new PhasedGenotype(2, 1));
            var genotype2 = Genotype.Create(new PhasedGenotype(2, 1));
            Assert.Equal(genotype1, genotype2);

            genotype2 = Genotype.Create(new PhasedGenotype(1, 2));
            Assert.NotEqual(genotype1, genotype2);

            genotype2 = Genotype.Create(new PhasedGenotype(1, 1));
            Assert.NotEqual(genotype1, genotype2);
        }
    }
}

