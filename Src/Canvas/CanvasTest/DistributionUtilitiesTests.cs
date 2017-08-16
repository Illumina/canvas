using System.Collections.Generic;
using System.Linq;
using CanvasCommon;
using Xunit;

namespace CanvasTest
{
    public class DistributionUtilitiesTests
    {
        [Fact]
        public void TestGetGenotypeCombinations()
        {
            int numSamples = 2;
            int altCN = 1;
            var result = DistributionUtilities.GetGenotypeCombinations(numSamples, altCN);

            Assert.Equal(new List<List<int>>
            {
                new []{ 1, 2}.ToList(),
                new []{ 2, 1}.ToList(),
            }, result);
        }
        
        [Fact]
        public void TestGetGenotypeCombinationsSingleSample()
        {
            int numSamples = 1;
            int altCN = 1;
            var result = DistributionUtilities.GetGenotypeCombinations(numSamples, altCN);

            Assert.Equal(new[]
            {
                new []{ 1 }.ToList(),
            }.ToList(), result);
        }

        [Fact]
        public void TestNegativeBinomialWrapperMaxDensity()
        {
            // maximal desnity should equal mean
            double mean = 50.0;
            double variance = 50.0;
            int maxValue = 200;
            int maxValueIndex = 49;
            var result = DistributionUtilities.NegativeBinomialWrapper(mean, variance, maxValue);
            Assert.Equal(maxValueIndex, result.FindIndex(x => x == result.Max()));
        }
    }
}