using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Principal;
using CanvasCommon;
using MathNet.Numerics.Distributions;
using System.Threading;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;

namespace CanvasPartition
{
    /// <summary>
    /// Multivariate distibutions for HMM
    /// </summary>
    public abstract class MultivariateDistribution
    {
        public abstract int GetDimensions();
        public abstract List<double> Mean();   // the Mean property
        public abstract double EstimateLikelihood(List<double> x);
    }

    public class MultivariatePoissonDistribution : MultivariateDistribution
    {
        private readonly List<Poisson> _poisson;

        public MultivariatePoissonDistribution(List<double> means)
        {
            _poisson = new List<Poisson>();
            foreach (double mean in means)
                _poisson.Add(new Poisson(mean));
        }

        public override int GetDimensions()
        {
            return _poisson.Count;
        }

        public override List<double> Mean()
        {
            return _poisson.Select(x => x.Lambda).ToList();
        }

        public void UpdatePoisson(double[] gamma, List<List<double>> data)
        {
            var m = this.GetDimensions();
            for (int dimension = 0; dimension < m; dimension++)
                _poisson[dimension] = new Poisson(CanvasCommon.Utilities.WeightedMean(data.Select(x => x[dimension]).ToList(), gamma.ToList()));
        }

        public override double EstimateLikelihood(List<double> x)
        {
            var m = this.GetDimensions();
            double likelihood = _poisson[0].Probability(Convert.ToInt32(x[0]));
            for (int i = 1; i < m; i++)
                likelihood *= _poisson[i].Probability(Convert.ToInt32(x[i]));
            if (Double.IsNaN(likelihood) || Double.IsInfinity(likelihood))
                likelihood = 0;
            return likelihood;
        }
    }

    public class MultivariateGaussianDistribution : MultivariateDistribution
    {
        private readonly Vector<double> _mean;
        private readonly Matrix<double> _covariance;

        public MultivariateGaussianDistribution(Vector<double> mean, Matrix<double> covariance)
        {
            _mean = mean;
            _covariance = covariance;
        }

        public override int GetDimensions()
        {
            return _mean.Count;
        }

        public override List<double> Mean()
        {
            return _mean.ToList();
        }
        public Matrix<double> Covariance   // the Mean property
        {
            get
            {
                return _covariance;
            }
        }
        public void UpdateMean(double[] gamma, List<List<double>> data)
        {
            var m = this.GetDimensions();
            for (int dimension = 0; dimension < m; dimension++)
                _mean[dimension] = CanvasCommon.Utilities.WeightedMean(data.Select(x => x[dimension]).ToList(), gamma.ToList());
        }

        /// <summary>
        /// Implements spherical covariance
        /// </summary>
        public void UpdateCovariance(double[] gamma, List<List<double>> data)
        {
            var m = this.GetDimensions();
            for (int dimension = 0; dimension < m; dimension++)
                _covariance[dimension, dimension] = CanvasCommon.Utilities.WeightedStandardDeviation(data.Select(x => x[dimension]).ToList(), gamma);
        }

        public override double EstimateLikelihood(List<double> x)
        {
            var m = this.GetDimensions();
            Vector<double> diff = Vector<double>.Build.Dense(m, 0.0);
            for (int i = 0; i < m; i++)
                diff[i] = x[i] - _mean[i];
            var exponent = -0.5 * diff.DotProduct(_covariance.Inverse() * diff);
            if (!Double.IsNaN(exponent)) //check for nans
            {
                var likelihood = 1.0 / (Math.Sqrt(2.0 * Math.PI * _covariance.Determinant())) * Math.Exp(exponent);
                if (Double.IsNaN(likelihood) || Double.IsInfinity(likelihood))
                    likelihood = 0;
                return likelihood;
            }
            return 0;
        }
    }

    /// <summary>
    /// Mixture of multivariate distibutions for HMM
    /// </summary>
    public abstract class MixtureDistibution
    {
        public abstract double EstimateLikelihood(List<double> x, int state);
        public abstract void WriteMeans();
        public abstract void UpdateMeans(double[][] gamma, List<List<double>> x);
        public abstract void UpdateCovariances(double[][] gamma, List<List<double>> x);
    }

    public class PoissonMixture : MixtureDistibution
    {

        private readonly List<MultivariatePoissonDistribution> _poissonDistributions;

        public override void UpdateCovariances(double[][] gamma, List<List<double>> x) {}

        public PoissonMixture(List<MultivariatePoissonDistribution> poissonDistributions)
        {
            _poissonDistributions = poissonDistributions;
        }

        public override void UpdateMeans(double[][] gamma, List<List<double>> x)
        {
            int stateCounter = 0;
            foreach (MultivariatePoissonDistribution poissonDistribution in _poissonDistributions)
            {
                poissonDistribution.UpdatePoisson(gamma[stateCounter], x);
                stateCounter++;
            }
        }

        public override double EstimateLikelihood(List<double> x, int state)
        {
            return _poissonDistributions[state].EstimateLikelihood(x);
        }


        public override void WriteMeans()
        {
            foreach (MultivariatePoissonDistribution gaussianDistribution in _poissonDistributions)
                for (int i = 0; i < gaussianDistribution.Mean().Count; i++)
                    Console.WriteLine($"Poisson mean for state {i} = {gaussianDistribution.Mean()[i]}");
        }
    }

    

    public class GaussianMixture : MixtureDistibution
    {

        private readonly List<MultivariateGaussianDistribution> _gaussianDistributions;
        public GaussianMixture(List<MultivariateGaussianDistribution> gaussianDistributions)
        {
            _gaussianDistributions = gaussianDistributions;
        }
        public override void UpdateMeans(double[][] gamma, List<List<double>> x)
        {
            int stateCounter = -1;
            foreach (MultivariateGaussianDistribution gaussianDistribution in _gaussianDistributions)
            {
                stateCounter++;
                if (stateCounter == 1)
                    continue;
                gaussianDistribution.UpdateMean(gamma[stateCounter], x);
            }
        }
        public override void WriteMeans()
        {
            foreach (MultivariateGaussianDistribution gaussianDistribution in _gaussianDistributions)
                for (int i = 0; i < gaussianDistribution.Mean().Count; i++)
                    Console.WriteLine($"Gaussian mean for state {i} = {gaussianDistribution.Mean()[i]}");
        }
        public override void UpdateCovariances(double[][] gamma, List<List<double>> x)
        {
            int stateCounter = -1;
            foreach (MultivariateGaussianDistribution gaussianDistribution in _gaussianDistributions)
            {
                stateCounter++;
                if (stateCounter == 1)
                    continue;
                gaussianDistribution.UpdateCovariance(gamma[stateCounter], x);
            }
        }
        public override double EstimateLikelihood(List<double> x, int state)
        {
            return _gaussianDistributions[state].EstimateLikelihood(x);
        }
        public int GetNumberOfStates()
        {
            return _gaussianDistributions.Count;
        }
    }
}
