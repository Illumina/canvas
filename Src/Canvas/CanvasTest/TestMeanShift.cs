using Xunit;
using Modified.Accord.MachineLearning;

namespace CanvasTest
{
    public class TestMeanShift
    {

        //This file/test exists to demonstrate the problem with the official Accord.Net MeanShift code.
        //If we get an update that fixes the issue, we can test it by changing DebugMeanShift
        //below to MeanShift.
        //Once it is working, this test can be eliminated.
        [Fact]
        public void TestClusterOrder()
        {
            // set of points with 4 obvious clusters
            double[][] points = {
                // plus,minus points
                new double[] { 11, -10},
                new double[] { 11, -12},
                new double[] { 10, -13},
                // plus,plus points
                new double[] { 10, 10},
                new double[] { 11, 13},
                new double[] { 10, 12},
                new double[] { 11, 10},
                // minus,plus points
                new double[] { -10, 10 },
                new double[] { -10, 11 },
                new double[] { -11, 10 },
                new double[] { -11, 11 },
                // minus,minus points
                new double[] { -10, -10 },
                new double[] { -11.5, -10 },
                new double[] { -13, -10 }
            };

            // for use in tests, count points in each cluster
            double minusMinusPointCount = 0;
            double minusPlusPointCount = 0;
            double plusMinusPointCount = 0;
            double plusPlusPointCount = 0;
            for (int i = 0; i < points.Length; ++i)
            {
                if (points[i][0] < 0 && points[i][1] < 0)
                {
                    minusMinusPointCount += 1;
                }
                else if (points[i][0] < 0 && points[i][1] > 0)
                {
                    minusPlusPointCount += 1;
                }
                else if (points[i][0] > 0 && points[i][1] < 0)
                {
                    plusMinusPointCount += 1;
                }
                else if (points[i][0] > 0 && points[i][1] > 0)
                {
                    plusPlusPointCount += 1;
                }
            }


            // MeanShift calculations
            Accord.Math.Random.Generator.Seed = 1;
            var meanShift = new MeanShift()
            {
                // Use a uniform kernel density
                Kernel = new Accord.Statistics.Distributions.DensityKernels.GaussianKernel(2),
                Bandwidth = 2
            };
            meanShift.UseParallelProcessing = false;

            var clustering = meanShift.Learn(points);

            int[] labels = clustering.Decide(points);


            // Test results

            // we should get 4 clusters
            Assert.True(clustering.Count == 4);

            // proportions of clusters should match, and point labels should assign them to modes
            // that make sense
            for (int i = 0; i < clustering.Count; ++i)
            {
                if (clustering.Modes[i][0] < 0 && clustering.Modes[i][1] < 0)
                {
                    Assert.Equal(clustering.Proportions[i], minusMinusPointCount / points.Length, 5);
                    for (int j = 0; j < points.Length; ++j)
                    {
                        if (labels[j] == i)
                            Assert.True(points[j][0] < 0 && points[j][1] < 0);
                    }
                }
                if (clustering.Modes[i][0] < 0 && clustering.Modes[i][1] > 0)
                {
                    Assert.Equal(clustering.Proportions[i], minusPlusPointCount / points.Length, 5);
                    for (int j = 0; j < points.Length; ++j)
                    {
                        if (labels[j] == i)
                            Assert.True(points[j][0] < 0 && points[j][1] > 0);
                    }
                }
                if (clustering.Modes[i][0] > 0 && clustering.Modes[i][1] < 0)
                {
                    Assert.Equal(clustering.Proportions[i], plusMinusPointCount / points.Length, 5);
                    for (int j = 0; j < points.Length; ++j)
                    {
                        if (labels[j] == i)
                            Assert.True(points[j][0] > 0 && points[j][1] < 0);
                    }
                }
                if (clustering.Modes[i][0] > 0 && clustering.Modes[i][1] > 0)
                {
                    Assert.Equal(clustering.Proportions[i], plusPlusPointCount / points.Length, 5);
                    for (int j = 0; j < points.Length; ++j)
                    {
                        if (labels[j] == i)
                            Assert.True(points[j][0] > 0 && points[j][1] > 0);
                    }
                }
            }
        }
    }
}