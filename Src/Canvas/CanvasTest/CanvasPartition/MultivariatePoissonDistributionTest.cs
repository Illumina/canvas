using Microsoft.VisualStudio.TestTools.UnitTesting;
using CanvasPartition;
using System.Collections.Generic;

namespace CanvasTest.CanvasPartition
{
    [TestClass()]
    public class MultivariatePoissonDistributionTestForZeroProbability
    {
        [TestMethod()]
        public void EstimateLikelihoodTest()
        {
            MultivariatePoissonDistribution test = new MultivariatePoissonDistribution(new List<double> { 1, 1 });
            double likelihood = test.EstimateLikelihood(new List<double> { 100000, 100000 });
            Assert.AreEqual(0, likelihood);
        }
    }
}