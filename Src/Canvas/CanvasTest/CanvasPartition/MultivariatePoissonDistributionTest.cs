using Xunit;
using CanvasPartition;
using System.Collections.Generic;

namespace CanvasTest.CanvasPartition
{
    public class MultivariatePoissonDistributionTestForZeroProbability
    {
        [Fact]
        public void EstimateLikelihoodTest()
        {
            MultivariatePoissonDistribution test = new MultivariatePoissonDistribution(new List<double> { 1, 1 });
            double likelihood = test.EstimateLikelihood(new List<double> { 100000, 100000 });
            Assert.Equal(0, likelihood);
        }
    }
}