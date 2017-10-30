using System;
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
            segmentOne.Maf = 0.5;

            var segmentTwo = new SegmentInfo();
            segmentTwo.Coverage = 1;
            segmentTwo.Maf = 0;

            var segments = new List<SegmentInfo>
            {
                segmentTwo,
                segmentTwo,
                segmentOne
            };
            var clustering = new DensityClusteringModel(segments, .5, 2, 2);
            clustering.EstimateDistance();
            double distanceThreshold = clustering.EstimateDc();
            clustering.GaussianLocalDensity(distanceThreshold);
            clustering.FindCentroids();
            var clusterCount = clustering.FindClusters(2);
            Assert.Equal(0, clusterCount);
        }

        [Fact]
        public void EstimateDc_tmpHigh_wrong()
        {
            List<double> distances = new List<double> { 6, 5, 4, 3, 2, 1 };
            double tmpLow = Double.MaxValue;
            double tmpHigh = Double.MinValue;
            /** SK: The logic is incorrect
                for Distance = {6,5,4,3,2,1}, tmpLow and tmpHigh would 1 and Double.MinValue, respectively
            **/

            foreach (double? element in distances)
            {
                if (element.HasValue && element < tmpLow && element > 0)
                    tmpLow = (double)element;
                else if (element.HasValue && element > tmpHigh)
                    tmpHigh = (double)element;
            }

            Assert.Equal(1D, tmpLow);
            Assert.Equal(tmpHigh, Double.MinValue);
        }
    }
}