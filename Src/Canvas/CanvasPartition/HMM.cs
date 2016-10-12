using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Principal;
using CanvasCommon;
using MathNet.Numerics.Distributions;
using System.Threading;

namespace CanvasPartition
{
    public class MultivariateGaussianDistribution
    {
        private readonly double[] _mean;
        private readonly double[][] _covariance;

        public MultivariateGaussianDistribution(double[] mean, double[][] covariance)
        {
            _mean = mean;
            _covariance = covariance;
        }

        public int GetDimensions()
        {
            return _mean.Length;
        }

        public double[] Mean   // the Mean property
        {
            get
            {
                return _mean;
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
                _covariance[dimension][dimension] = CanvasCommon.Utilities.WeightedStandardDeviation(data.Select(x => x[dimension]).ToList(), gamma);
        }

        public double EstimateLikelihood(List<double> x)
        {
            var m = this.GetDimensions();
            double[] diff = new double[m];
            for (int i = 0; i < m; i++)
                diff[i] = x[i] - _mean[i];
            var inverseCovariance = MatrixUtilities.MatrixInverse(_covariance);
            var tempval = MatrixUtilities.MatrixVectorProduct(inverseCovariance, diff);
            var exponent = -0.5 * CanvasCommon.Utilities.DotProduct(diff, tempval);
            if (!Double.IsNaN(exponent)) //check for nans
            {
                var likelihood = 1.0 / (Math.Sqrt(2.0 * Math.PI * MatrixUtilities.MatrixDeterminant(_covariance))) *
                                 Math.Exp(exponent);
                if (Double.IsNaN(likelihood))
                    likelihood = 0;
                return likelihood;
            }
            return 0;
        }
    }

    public class GaussianMixture
    {

        private readonly List<MultivariateGaussianDistribution> _gaussianDistributions;
        public GaussianMixture(List<MultivariateGaussianDistribution> gaussianDistributions)
        {
            _gaussianDistributions = gaussianDistributions;
        }
        public void UpdateMeans(double[][] gamma, List<List<double>> x)
        {
            int stateCounter = 0;
            foreach (MultivariateGaussianDistribution gaussianDistribution in _gaussianDistributions)
            {
                gaussianDistribution.UpdateMean(gamma[stateCounter], x);
                stateCounter++;
            }

        }
        public void UpdateCovariances(double[][] gamma, List<List<double>> x)
        {
            int stateCounter = 0;
            foreach (MultivariateGaussianDistribution gaussianDistribution in _gaussianDistributions)
            {
                gaussianDistribution.UpdateCovariance(gamma[stateCounter], x);
                stateCounter++;
            }
        }
        public double EstimateLikelihood(List<double> x, int state)
        {
            return _gaussianDistributions[state].EstimateLikelihood(x);
        }
        public int GetNumberOfStates()
        {
            return _gaussianDistributions.Count;
        }
    }

    public class HiddenMarkovModel
    {
        // hidden variable parameters
        private readonly GaussianMixture _emission;
        private readonly double[][] _transition;
        private readonly double[] _stateProbabilities;
        // alpha-beta algorithm
        private readonly double[][] _beta;
        private readonly double[][] _alpha;
        // marginal posterior distibution
        private readonly double[][] _gamma;
        // joint posterior distibution
        private readonly double[][][] _epsilon;
        // dimention variables
        public int nStates;
        public int T;
        public const double selfTransition = 0.9;


        public HiddenMarkovModel(List<List<double>> data, List<MultivariateGaussianDistribution> gaussianMixtures)
        {
            // HMM set-up
            nStates = gaussianMixtures.Count;
            T = data.Count;
            _stateProbabilities = new double[nStates];
            _emission = new GaussianMixture(gaussianMixtures);
            _transition = MatrixUtilities.MatrixCreate(nStates, nStates);
            for (int i = 0; i < nStates; i++) { 
                for (int j = 0; j < nStates; j++)
                {
                    if (i == j)
                        _transition[i][j] = selfTransition;
                    else
                        _transition[i][j] = (1.0 - selfTransition)/(nStates - 1);
                }
            }

            // alpha-beta algorithm
            _beta = MatrixUtilities.MatrixCreate(T, nStates);
            _alpha = MatrixUtilities.MatrixCreate(T, nStates);
            // marginal posterior distibution
            _gamma = MatrixUtilities.MatrixCreate(nStates, T);
            // joint posterior distibution
            _epsilon = MatrixUtilities.MatrixCreate(T, nStates, nStates);
        }

        public void IncreaseVariance(int state = 2, double multiplicationFactor = 2.0)
        {
            _transition[state][state] = _transition[state][state]*multiplicationFactor;
        }

        public void UpdateTransition()
        {
            double [][] numerator = MatrixUtilities.MatrixCreate(nStates, nStates);
            double [] denominator = new double[nStates];
            for (int t = 0; t < T-1; t++)
            {
                for (int i = 0; i < nStates; i++)
                {
                    for (int j = 0; j < nStates; j++)
                        numerator[i][j] += _epsilon[t][i][j];
                    denominator[i] += _gamma[i][t];
                }
            }
            // update transition and initial state probabilities
            for (int i = 0; i < nStates; i++)
                denominator[i] = Math.Max(denominator[i], Single.MinValue);           
            
            for (int i = 0; i < nStates; i++)
            {
                for (int j = 0; j < nStates; j++)
                    _transition[i][j] = numerator[i][j] / denominator[i];
                _stateProbabilities[i] = denominator[i] / denominator.Average();
            }

        }

        public void UpdateEmission(List<List<double>> x)
        {
            _emission.UpdateCovariances(_gamma, x);
            _emission.UpdateMeans(_gamma, x);
        }

        public double Forward(List<List<double>> x)
        {
            // Initialization 
            double likelihood = 0;
            for (int j = 0; j < nStates; j++)
                _alpha[0][j] = 1/(double)nStates;

            // Induction 
            for (int t = 0; t < T-1; t++)
            {
                double scaler = 0;
                for (int j = 0; j < nStates; j++)
                {
                    for (int i = 0; i < nStates; i++)
                        _alpha[t + 1][j] += _alpha[t][i] * _transition[i][j];

                    _alpha[t + 1][j] *= _emission.EstimateLikelihood(x[t + 1], j);
                    scaler += _alpha[t + 1][j];
                }
                for (int j = 0; j < nStates; j++)
                {
                    _alpha[t + 1][j] = _alpha[t + 1][j] / scaler;
                }
                likelihood += scaler;
            }
            return -Math.Log(likelihood);
        }
        public void Backward(List<List<double>> x)
        {
            // Initialization 
            for (int j = 0; j < nStates; j++)
                _beta[T-1][j] = 1 / (double)nStates;


            // Induction 
            int t = T - 2;
            while (t >= 0)
            {
                double scaler = 0;
                for (int i = 0; i < nStates; i++)
                {
                    for (int j = 0; j < nStates; j++)
                        _beta[t][i] += _transition[i][j] * _emission.EstimateLikelihood(x[t + 1], j) * _beta[t + 1][j];
                    scaler += _beta[t][i];
                }
                for (int i = 0; i < nStates; i++)
                    _beta[t][i] = _beta[t][i] / scaler;
                t--;
            }
        }


        public void MaximisationStep(List<List<double>> x)
        {
            UpdateTransition();
            UpdateEmission(x);
        }

        public double EstimationStep(List<List<double>> x)
        {
            List<ThreadStart> tasks = new List<ThreadStart>();
            double likelihood = 0;
            tasks.Add(new ThreadStart(() => { likelihood = Forward(x); }));
            tasks.Add(new ThreadStart(() => { Backward(x); }));
            Isas.Shared.Utilities.Utilities.DoWorkParallelThreads(tasks);

            for (int t = 0; t < T-1; t++)
            {
                for (int i = 0; i < nStates; i++)
                {
                    for (int j = 0; j < nStates; j++) { 
                        _epsilon[t][i][j] = _alpha[t][i] * _transition[i][j] * _emission.EstimateLikelihood(x[t + 1], j) *
                            _beta[t + 1][j] / likelihood;
                        _gamma[i][t] += _epsilon[t][i][j];
                    }
                }
            }
            return likelihood;
        }

        public void FindMaximalLikelyhood(List<List<double>> x, string chr)
        {
            double likelihoodDifferenceThreshold = 0.01;
            const int maxIterations = 5;
            List<double> likelihoods = new List<double>();
            double oldLikelihood = EstimationStep(x);
            likelihoods.Add(oldLikelihood);
            MaximisationStep(x);
            int iteration = 1;
            double likelihoodDifference = Double.MaxValue;
            while (iteration < maxIterations && likelihoodDifference > likelihoodDifferenceThreshold)
            {
                Console.WriteLine($"{DateTime.Now} EM step {iteration} for chromosome {chr}");
                var newLikelihood = EstimationStep(x);
                MaximisationStep(x);
                likelihoodDifference = Math.Abs(oldLikelihood - newLikelihood);
                likelihoods.Add(newLikelihood);
                oldLikelihood = newLikelihood;
                iteration++;
            }
            Console.WriteLine($"Completed EM step for chromosome {chr}");

        }

        public List<int> BestPathViterbi(List<List<double>> x)
        {
            // Initialization 
            double[][] bestScore = MatrixUtilities.MatrixCreate(T + 1, nStates);
            int[][] bestStateSequence = new int[T + 1][];
            for (int i = 0; i < T; ++i)
                bestStateSequence[i] = new int[nStates];

            for (int j = 0; j < nStates; j++) { 
                bestScore[0][j] = this._stateProbabilities[j];
                bestStateSequence[0][j] = 0;
            }

            // Induction 
            for (int t = 1; t < T - 1; t++)
            {
                for (int j = 0; j < nStates; j++)
                {
                    int state = 0;
                    double max = Double.MinValue;
                    for (int i = 0; i < nStates; i++)
                    {
                        var tmpMax = bestScore[t-1][i] + Math.Log(_transition[i][j]) + Math.Log(_emission.EstimateLikelihood(x[t], j));
                        if (tmpMax > max)
                        {
                            state = i;
                            max = tmpMax;
                        }
                    }
                    bestScore[t][j] = max;
                    bestStateSequence[t][j] = state;
                }
            }

            var backtrack = T-1;
            int bestState = 0;
            List<int> bestStates = new List<int>(T+1);
            for (int i = 0; i < nStates; i++)
            {
                double max = Double.MinValue;
                var tmpMax = bestScore[T][i];
                if (tmpMax > max)
                {
                    bestState = i;
                    max = tmpMax;
                }
            }
            bestStates.Add(backtrack);
            while (backtrack >= 0)
            {
                bestStates.Add(bestStateSequence[backtrack][bestState]);
                backtrack--;
            }
            bestStates.Reverse();
            return bestStates;
        }
    }
}
