using System;
using System.Collections.Generic;
using System.Linq;
using CanvasCommon;
using MathNet.Numerics.Distributions;
using Combinatorics.Collections;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics;

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

    public class MultivariateNegativeBinomial : MultivariateDistribution
    {
        private readonly List<List<double>> _negativeBinomials;
        public List<double> Means { get; set; }
        public List<double> Variances { get; set; }


        public MultivariateNegativeBinomial(List<double> means, List<double> variances, int maxValue)
        {
            _negativeBinomials = new List<List<double>>();
            for (int i = 0; i < means.Count; i++)
                _negativeBinomials.Add(DistibutionUtilities.NegativeBinomialWrapper(means[i], variances[i], maxValue));
            Means = means;
            Variances = variances;
        }

        public List<List<double>> NegativeBinomials => _negativeBinomials;

        public override int GetDimensions()
        {
            return NegativeBinomials.Count;
        }

        public override List<double> Mean()
        {
            return Means;
        }

        public void UpdateNegativeBinomial(double[] gamma, List<List<double>> data, List<double> variance, int maxValue)
        {
            var m = this.GetDimensions();
            for (int dimension = 0; dimension < m; dimension++)
            {
                // var newVariance = CanvasCommon.Utilities.WeightedVariance(data.Select(x => x[dimension]).ToList(), gamma);
                var newMeans = CanvasCommon.Utilities.WeightedMean(data.Select(x => x[dimension]).ToList(),
                    gamma.ToList());
                Means[dimension] = newMeans;
                NegativeBinomials[dimension] = DistibutionUtilities.NegativeBinomialWrapper(newMeans, variance[dimension], maxValue);
            }
        }
                
        public override double EstimateLikelihood(List<double> x)
        {
            var m = this.GetDimensions();
            double likelihood = NegativeBinomials[0][Convert.ToInt32(x[0])];
            for (int i = 1; i < m; i++)
                likelihood *= NegativeBinomials[i][Convert.ToInt32(x[i])];
            if (Double.IsNaN(likelihood) || Double.IsInfinity(likelihood))
                likelihood = 0;
            return likelihood;
        }

        public double Probability(int cnState, int pointCoverage)
        {
            return NegativeBinomials[cnState][pointCoverage];
        }
    }

    public class MultivariatePoissonDistribution : MultivariateDistribution
    {
        public MultivariatePoissonDistribution(List<double> means)
        {
            Poisson = new List<Poisson>();
            foreach (double mean in means)
                Poisson.Add(new Poisson(mean));
        }

        public List<Poisson> Poisson { get; }

        public override int GetDimensions() => Poisson.Count;

        public override List<double> Mean()
        {
            return Poisson.Select(x => x.Lambda).ToList();
        }

        public void UpdatePoisson(double[] gamma, List<List<double>> data)
        {
            var m = this.GetDimensions();
            for (int dimension = 0; dimension < m; dimension++)
                Poisson[dimension] = new Poisson(CanvasCommon.Utilities.WeightedMean(data.Select(x => x[dimension]).ToList(), gamma.ToList()));
        }

        public override double EstimateLikelihood(List<double> x)
        {
            var m = this.GetDimensions();
            double likelihood = Poisson[0].Probability(Convert.ToInt32(x[0]));
            for (int i = 1; i < m; i++)
                likelihood *= Poisson[i].Probability(Convert.ToInt32(x[i]));
            if (Double.IsNaN(likelihood) || Double.IsInfinity(likelihood))
                likelihood = 0;
            return likelihood;
        }
    }

    public class MultivariateGaussianDistribution : MultivariateDistribution
    {
        private readonly Vector<double> _mean;

        public MultivariateGaussianDistribution(Vector<double> mean, Matrix<double> covariance)
        {
            _mean = mean;
            Covariance = covariance;
        }

        public override int GetDimensions()
        {
            return _mean.Count;
        }

        public override List<double> Mean()
        {
            return _mean.ToList();
        }
        public Matrix<double> Covariance { get; }

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
                Covariance[dimension, dimension] = CanvasCommon.Utilities.WeightedStandardDeviation(data.Select(x => x[dimension]).ToList(), gamma);
        }

        public override double EstimateLikelihood(List<double> x)
        {
            var m = this.GetDimensions();
            Vector<double> diff = Vector<double>.Build.Dense(m, 0.0);
            for (int i = 0; i < m; i++)
                diff[i] = x[i] - _mean[i];
            var exponent = -0.5 * diff.DotProduct(Covariance.Inverse() * diff);
            if (!Double.IsNaN(exponent)) //check for nans
            {
                var likelihood = 1.0 / (Math.Sqrt(2.0 * Math.PI * Covariance.Determinant())) * Math.Exp(exponent);
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
        public abstract void UpdateMeans(double[][] gamma, List<List<double>> x, List<List<double>> variance = null);
        public abstract void UpdateCovariances(double[][] gamma, List<List<double>> x);
        public abstract double EstimateViterbiLikelihood(List<double> x, int currentCnState, List<double> haploidMeans, double[] transition);

    }

    public static class Utilities
    {
        public static List<List<int>> GetGenotypeCombinations(int numberOfDimensions, int currentState)
        {
            if (numberOfDimensions <= 0)
                throw new ArgumentOutOfRangeException(nameof(numberOfDimensions));
            const int diploidState = 2;
            const int numbnerOfStates = 2;
            var upperSetBound = Math.Pow(numbnerOfStates, numberOfDimensions);
            var allCombinations = new List<List<int>>(Convert.ToInt32(upperSetBound));
            for (int numberOfDiploidStates = 1; numberOfDiploidStates < numberOfDimensions; numberOfDiploidStates++)
            {
                var states = Enumerable.Repeat(currentState, numberOfDimensions - numberOfDiploidStates)
                    .Concat(Enumerable.Repeat(diploidState, numberOfDiploidStates));
                var permutations = new Permutations<int>(states.ToList(), GenerateOption.WithoutRepetition);
                var list = permutations.Select(x => x.ToList()).ToList();
                allCombinations.AddRange(list);
            }
            return allCombinations;
        }

        public static List<double> NegativeBinomialWrapper(double mean, double variance, int maxValue)
        {
            var density = Enumerable.Repeat(0.0, maxValue).ToList();
            double r = Math.Pow(Math.Max(mean, 0.1), 2) / (Math.Max(variance, mean*1.2) - mean);
            for (int x = 0; x < maxValue; x++)
            {
                var tmpDensity = Math.Exp(Math.Log(Math.Pow(1 + mean / r, -r)) + Math.Log(Math.Pow(mean / (mean + r), x)) +  SpecialFunctions.GammaLn(r + x) - 
                             SpecialFunctions.FactorialLn(x) - SpecialFunctions.GammaLn(r));
                density[x] = Double.IsNaN(tmpDensity) || Double.IsInfinity(tmpDensity) ? 0 : tmpDensity;
            }              
            return density;
        }
    }

    public class NegativeBinomialMixture : MixtureDistibution
    {
        private readonly List<MultivariateNegativeBinomial> _negativeBinomialDistributions;
        public Dictionary<int, List<List<int>>> GenotypePermutations { get; }
        public override void UpdateCovariances(double[][] gamma, List<List<double>> x) { }

        public NegativeBinomialMixture(List<MultivariateNegativeBinomial> negativeBinomialDistributions, List<double> haploidMeans)
        {
            _negativeBinomialDistributions = negativeBinomialDistributions;
            GenotypePermutations = new Dictionary<int, List<List<int>>>();

            for (int cnState = 0; cnState < negativeBinomialDistributions.Count; cnState++)
            {
                List<List<int>> currentGenotypePermutation = null;
                int nDimensions = _negativeBinomialDistributions[cnState].GetDimensions();
                if (_negativeBinomialDistributions[cnState].Mean().Average() < haploidMeans.Average() * 1.5 ||
                    _negativeBinomialDistributions[cnState].Mean().Average() > haploidMeans.Average() * 1.5)
                {
                    currentGenotypePermutation = DistibutionUtilities.GetGenotypeCombinations(nDimensions, cnState);
                }
                GenotypePermutations[cnState] = currentGenotypePermutation;
            }
        }

        public override void UpdateMeans(double[][] gamma, List<List<double>> x, List<List<double>> variance)
        {
            int stateCounter = 0;
            var maxValue = x.Select(y => Convert.ToInt32(y.Max())).Max() + 10;
            foreach (MultivariateNegativeBinomial poissonDistribution in _negativeBinomialDistributions)
            {
                poissonDistribution.UpdateNegativeBinomial(gamma[stateCounter], x, variance[stateCounter], maxValue);
                stateCounter++;
            }
        }

        public override double EstimateLikelihood(List<double> x, int state)
        {
            return _negativeBinomialDistributions[state].EstimateLikelihood(x);
        }

        public override double EstimateViterbiLikelihood(List<double> data, int currentCnState, List<double> haploidMeans, double[] transition)
        {
            if (GenotypePermutations[currentCnState] == null)
                return Math.Log(_negativeBinomialDistributions[currentCnState].EstimateLikelihood(data)) + Math.Log(transition[currentCnState]);

            double maxLikelyhood = Double.MinValue;
            var bestState = new List<int>();

            foreach (List<int> currentGenotypePermutations in GenotypePermutations[currentCnState])
            {
                double emissionLikelihood = 1.0;
                int cnState = 0;
                foreach (int cnGenotype in currentGenotypePermutations)
                {
                    var pointCoverage = Convert.ToInt32(data[cnState]);
                    if (cnGenotype == 0 | cnGenotype == 1)
                        emissionLikelihood *= Math.Max(_negativeBinomialDistributions[0].Probability(cnState, pointCoverage),
                                _negativeBinomialDistributions[1].Probability(cnState, pointCoverage));
                    else if (cnGenotype == 3 | cnGenotype == 4)
                        emissionLikelihood *= Math.Max(_negativeBinomialDistributions[3].Probability(cnState, pointCoverage),
                                _negativeBinomialDistributions[4].Probability(cnState, pointCoverage));
                    else
                        emissionLikelihood *= _negativeBinomialDistributions[cnGenotype].Probability(cnState, pointCoverage);
                    cnState++;
                }
                if (Double.IsNaN(emissionLikelihood) || Double.IsInfinity(emissionLikelihood))
                    return 0;
                if (maxLikelyhood < emissionLikelihood)
                {
                    bestState = currentGenotypePermutations;
                    maxLikelyhood = emissionLikelihood;

                }
            }
            return Math.Log(maxLikelyhood) + Math.Log(bestState.Select(state => transition[state]).Average());
        }


        public override void WriteMeans()
        {
            foreach (MultivariateNegativeBinomial negativeBinomial in _negativeBinomialDistributions)
                for (int i = 0; i < negativeBinomial.Mean().Count; i++) 
                    Console.WriteLine($"Negative Binomial mean/variance for state {i} = {negativeBinomial.Means[i]} : {negativeBinomial.Variances[i]}");
        }
    }

    public class PoissonMixture : MixtureDistibution
    {

        private readonly List<MultivariatePoissonDistribution> _poissonDistributions;
        public Dictionary<int, List<List<int>>> GenotypePermutations { get; }


        public override void UpdateCovariances(double[][] gamma, List<List<double>> x) { }

        public PoissonMixture(List<MultivariatePoissonDistribution> poissonDistributions, List<double> haploidMeans)
        {
            _poissonDistributions = poissonDistributions;
            GenotypePermutations = new Dictionary<int, List<List<int>>>();

            for (int cnState = 0; cnState < _poissonDistributions.Count; cnState++)
            {
                List<List<int>> currentGenotypePermutation = null;
                int nDimensions = _poissonDistributions[cnState].GetDimensions();
                if (_poissonDistributions[cnState].Mean().Average() < haploidMeans.Average() * 1.5 ||
                    _poissonDistributions[cnState].Mean().Average() > haploidMeans.Average() * 1.5)
                {
                    currentGenotypePermutation = DistibutionUtilities.GetGenotypeCombinations(nDimensions, cnState);
                }
                GenotypePermutations[cnState] = currentGenotypePermutation;
            }
        }

        public override void UpdateMeans(double[][] gamma, List<List<double>> x, List<List<double>> variance = null)
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

        public override double EstimateViterbiLikelihood(List<double> x, int currentCnState, List<double> haploidMeans, double[] transition)
        {
            if (GenotypePermutations[currentCnState] != null)
            {
                double maxLikelyhood = Double.MinValue;
                var bestState = new List<int>();

                foreach (List<int> states in GenotypePermutations[currentCnState])
                {
                    double emissionLikelihood = 1.0;
                    int stateCounter = 0;
                    foreach (int state in states)
                    {
                        emissionLikelihood *= _poissonDistributions[state].Poisson[stateCounter].Probability(Convert.ToInt32(x[stateCounter]));
                        stateCounter++;
                    }
                    if (Double.IsNaN(emissionLikelihood) || Double.IsInfinity(emissionLikelihood))
                        emissionLikelihood = 0;
                    if (maxLikelyhood < emissionLikelihood)
                    {
                        bestState = states;
                        maxLikelyhood = emissionLikelihood;

                    }
                }
                return Math.Log(maxLikelyhood) + Math.Log(bestState.Select(state => transition[state]).Average());
            }

            return Math.Log(_poissonDistributions[currentCnState].EstimateLikelihood(x)) + Math.Log(transition[currentCnState]);
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
        public override void UpdateMeans(double[][] gamma, List<List<double>> x, List<List<double>> variance = null)
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

        public override double EstimateViterbiLikelihood(List<double> x, int currentCnState, List<double> haploidMeans, double[] transition)
        {
            throw new NotImplementedException();
        }
    }
}
