using System.Collections.Generic;
using CanvasCommon;
using Xunit;

namespace CanvasTest
{
    public class DensityClusteringTests
    {
        [Fact]
        public void FindClusters_NoSegments_NoClusters()
        {
            var segmentOne = new SegmentInfo();
            segmentOne.Coverage = 2;
            segmentOne.MAF = 0.5;

            var segmentTwo = new SegmentInfo();
            segmentTwo.Coverage = 1;
            segmentTwo.MAF = 0;

            var segments = new List<SegmentInfo>
            {
                segmentOne,
                segmentTwo,
                segmentTwo
            };
            var clustering = new DensityClusteringModel(segments, .5, 2, 2);
            clustering.EstimateDistance();
            double distanceThreshold = clustering.EstimateDc();
            clustering.GaussianLocalDensity(distanceThreshold);
            clustering.FindCentroids();
            var clusterCount = clustering.FindClusters(2);
            Assert.Equal(2, clusterCount);
        }
    }
}