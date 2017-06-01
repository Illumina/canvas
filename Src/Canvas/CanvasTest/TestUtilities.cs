using System;
using System.Collections.Generic;
using System.Linq;
using CanvasCommon;
using Illumina.Common;
using Xunit;

namespace CanvasTest
{
    public class TestUtilities
    {
        [Theory]
        [InlineData("weightedaverage", CanvasNormalizeMode.WeightedAverage, false)]
        [InlineData("WeightedAverage", CanvasNormalizeMode.WeightedAverage, false)]
        [InlineData("bestlr2", CanvasNormalizeMode.BestLR2, false)]
        [InlineData("BestLR2", CanvasNormalizeMode.BestLR2, false)]
        [InlineData("pca", CanvasNormalizeMode.PCA, false)]
        [InlineData("PCA", CanvasNormalizeMode.PCA, false)]
        [InlineData("badmode", null, true)]
        public void TestParseCanvasNormalizeMode(string modeStr, CanvasNormalizeMode expectedMode, bool expectException)
        {
            if (expectException)
            {
                Exception ex = Assert.Throws<Exception>(() => Utilities.ParseCanvasNormalizeMode(modeStr));
                Assert.Equal(String.Format("Invalid CanvasNormalize mode '{0}'", modeStr), ex.Message);
            }
            else
            {
                CanvasNormalizeMode mode = Utilities.ParseCanvasNormalizeMode(modeStr);
                Assert.Equal(expectedMode, mode);
            }
        }

        [Theory]
        [InlineData(-5, 5)]
        [InlineData(0, 5)]
        [InlineData(-5, 0)]
        public void TestGoldenSectionSearch(double a, double b)
        {
            double res = Utilities.GoldenSectionSearch(x => x * x, a, b);
            Assert.True(Math.Abs(res) < 0.001);
        }

        [Fact]
        public void TestIsSubset()
        {
            string[] a = new string[] { "1", "2" };
            string[] b = new string[] { "1", "2", "3" };

            Assert.True(Utilities.IsSubset(a, b));
            Assert.False(Utilities.IsSubset(b, a));
        }

        [Fact]
        public void TestTwoNorm()
        {
            double[] v = new double[10];
            for (int i = 0; i < v.Length; i++) { v[i] = 1; }

            double size = Utilities.TwoNorm(v);
            double expectedSize = Math.Sqrt(10);

            Assert.True(Math.Abs(size - expectedSize) < 0.001);
        }

        [Theory]
        [InlineData(0, 0)]
        [InlineData(1, 1)]
        public void TestNormalizeBy2Norm(double value, double expected)
        {
            double[] v = new double[10];
            for (int i = 0; i < v.Length; i++) { v[i] = value; }

            double[] u = Utilities.NormalizeBy2Norm(v);
            for (int i = 0; i < u.Length; i++)
            {
                Assert.True(Math.Abs(u[i] - expected / Math.Sqrt(v.Length)) < 0.001);
            }
        }

        [Fact]
        public void TestDotProduct()
        {
            double[] v = new double[10];
            for (int i = 0; i < v.Length; i++) { v[i] = 1; }

            double dotProduct = Utilities.DotProduct(v, v);

            Assert.True(Math.Abs(dotProduct - v.Length) < 0.001);

            double[] w = new double[5];
            for (int i = 0; i < w.Length; i++) { w[i] = 1; }
            IlluminaException ex = Assert.Throws<IlluminaException>(() => Utilities.DotProduct(v, w));
            Assert.Equal("Vectors must be of the same dimension to calculate dot product.", ex.Message);
        }

        [Fact]
        public void TestAreOrthogonal()
        {
            double[] v1 = new double[] { 1, 1 };
            double[] v2 = new double[] { 1, -1 };

            Assert.False(Utilities.AreOrthogonal(v1, v1));
            Assert.True(Utilities.AreOrthogonal(v1, v2));
        }

        [Fact]
        public void TestProject()
        {
            double[] axis = new double[10];
            for (int i = 0; i < axis.Length; i++) { axis[i] = 1 / Math.Sqrt(axis.Length); }
            double[] u = new double[10];
            u[0] = 1;

            double[] w = Utilities.Project(u, axis);
            for (int i = 0; i < w.Length; i++)
            {
                Assert.True(Math.Abs(w[i] - 1.0 / axis.Length) < 0.001);
            }

            IlluminaException ex = Assert.Throws<IlluminaException>(() => Utilities.Project(u, new double[5]));
            Assert.Equal("Vector and the axis must be of the same dimension.", ex.Message);
        }

        public static List<double[]> GetHelmertBasis(int p)
        {
            List<double[]> axes = new List<double[]>();
            // Helmert transformation
            for (int i = 1; i <= p; i++)
            {
                double[] v = new double[p];
                double size = Math.Sqrt(i < p ? (i + i * i) : i);
                for (int j = 0; j < i; j++)
                {
                    v[j] = 1 / size;
                }
                if (i < p)
                    v[i] = -i / size;
                axes.Add(v);
            }

            return axes;
        }

        [Fact]
        public void TestProject2()
        {
            List<double[]> axes = GetHelmertBasis(10);

            // project u
            double[] u = new double[10];
            for (int i = 0; i < u.Length; i++)
            {
                u[i] = i;
            }

            double[] w = Utilities.Project(u, axes); // w should be the same as u after projection
            for (int i = 0; i < w.Length; i++)
            {
                Assert.True(Math.Abs(u[i] - w[i]) < 0.001);
            }

            IlluminaException ex = Assert.Throws<IlluminaException>(() => Utilities.Project(u, (IEnumerable<double[]>)null));
            Assert.Equal("No axes to project onto.", ex.Message);
            ex = Assert.Throws<IlluminaException>(() => Utilities.Project(u, new List<double[]>()));
            Assert.Equal("No axes to project onto.", ex.Message);

            ex = Assert.Throws<IlluminaException>(() => Utilities.Project(u, new List<double[]>() { new double[5] }));
            Assert.Equal("Vector and the axes must be of the same dimension.", ex.Message);
        }

        [Fact]
        public void TestProjectToComplement()
        {
            List<double[]> axes = GetHelmertBasis(10);

            // project u
            double[] u = new double[10];
            for (int i = 0; i < u.Length; i++)
            {
                u[i] = i;
            }

            double[] w = Utilities.ProjectToComplement(u, axes); // w should be 0 after projection
            for (int i = 0; i < w.Length; i++)
            {
                Assert.True(Math.Abs(w[i]) < 0.001);
            }

            double[] x = Utilities.ProjectToComplement(u, axes.Take(5));
            double[] y = Utilities.Project(u, axes.Take(5));
            double dotProduct = Utilities.DotProduct(x, y); // should be 0
            Assert.True(Math.Abs(dotProduct) < 0.001);
        }

        [Fact]
        public void TestMedianFilter()
        {
            float[] values = new float[] { 2, 1, 3, 5, 4, 6, 7, 8 };
            float[] expected = new float[] { 1.5f, 2, 3, 4, 5, 6, 7, 7.5f };
            var smoothedValues = Utilities.MedianFilter(values, 1).ToArray();

            for (int i = 0; i < values.Length; i++)
            {
                Assert.Equal(expected[i], smoothedValues[i]);
            }
        }
    }
}