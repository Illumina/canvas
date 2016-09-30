using System;
using System.Collections.Generic;
using System.Linq;

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


        public void UpdateMean(double[] gamma, List<List<double>> x)
        {
            var n = gamma.Length;
            var m = this.GetDimensions();
            double[] numerator = new double[n];
            double denominator = 0;

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; i < m; i++)
                    numerator[i] += x[i][j] * gamma[i];
                denominator += gamma[i];
            }
            // update mean
            for (int i = 0; i < gamma.Length; i++)
                _mean[i] = numerator[i] / denominator;
        }

        /// <summary>
        /// Implements spherical covariance
        /// </summary>
        public void UpdateCovariance(double[] gamma, List<List<double>> x)
        {
            var n = gamma.Length;
            var m = this.GetDimensions();
            double[] numerator = new double[n];
            double denominator = 0;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; i < m; i++)
                    numerator[i] += Math.Pow(x[i][j] - _mean[j], 2) * gamma[i];
                denominator += gamma[i];
            }
            for (int i = 0; i < gamma.Length; i++)
                _covariance[i][i] = numerator[i] / denominator;
        }

        public double EstimateLikelihood(List<double> x)
        {
            var m = this.GetDimensions();
            double[] diff = new double[m];
            for (int i = 0; i < m; i++)
                diff[i] = x[i] - _mean[i];
            var inverseCovariance = CanvasCommon.Utilities.MatrixInverse(_covariance);
            var tempval = CanvasCommon.Utilities.MatrixVectorProduct(inverseCovariance, diff);
            var exponent = CanvasCommon.Utilities.DotProduct(diff, tempval);
            if (!Double.IsNaN(exponent)) //check for nans
            {
                var likelihood = 1.0 / (Math.Sqrt(2.0 * Math.PI * CanvasCommon.Utilities.MatrixDeterminant(_covariance))) *
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


        public HiddenMarkovModel(int nStates, int T)
        {
            // HMM set-up
            _stateProbabilities = new double[nStates];
            _emission = new GaussianMixture(new List<MultivariateGaussianDistribution>());
            _transition = CanvasCommon.Utilities.MatrixCreate(nStates, nStates);
            // alpha-beta algorithm
            _beta = CanvasCommon.Utilities.MatrixCreate(T, nStates);
            _alpha = CanvasCommon.Utilities.MatrixCreate(T, nStates);
            // marginal posterior distibution
            _gamma = CanvasCommon.Utilities.MatrixCreate(T, nStates);
            // joint posterior distibution
            _epsilon = CanvasCommon.Utilities.MatrixCreate(T, nStates, nStates);
            this.nStates = nStates;
            this.T = T;
        }

        public void UpdateTransition()
        {
            double [][] numerator = CanvasCommon.Utilities.MatrixCreate(nStates, nStates);
            double [] denominator = new double[nStates];
            for (int t = 0; t < T; t++)
            {
                for (int i = 0; i < nStates; i++)
                {
                    for (int j = 0; j < nStates; i++)
                        numerator[i][j] += _epsilon[t][i][j];
                    denominator[i] = _gamma[t][i];
                }
            }
            // update transition and initial state probabilities
            for (int i = 0; i < nStates; i++)
            {
                for (int j = 0; j < nStates; i++)
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
            double scaler = 0;
            double likelihood = 0;
            for (int j = 0; j < nStates; j++)
            {
                _alpha[0][j] *= _emission.EstimateLikelihood(x[0], j);
                scaler += _alpha[0][j];
            }
            for (int j = 0; j < nStates; j++)
                _alpha[0][j] = _alpha[0][j] / scaler;

            // Induction 
            for (int t = 1; t < T; t++)
            {
                scaler = 0;
                for (int j = 0; j < nStates; j++)
                {
                    for (int i = 0; i < nStates; i++)
                        _alpha[t + 1][j] = _alpha[t][i] * _transition[i][j];

                    _alpha[t + 1][j] *= _emission.EstimateLikelihood(x[t + 1], j);
                    scaler += _alpha[t + 1][j];
                }
                for (int j = 0; j < nStates; j++)
                {
                    _alpha[t + 1][j] = _alpha[t + 1][j] / scaler;
                }
                likelihood += scaler;
            }
            return likelihood;
        }
        public void Backward(List<List<double>> x)
        {
            // Initialization 
            for (int j = 0; j < nStates; j++)
                _beta[0][j] = 1 / (double)nStates;


            // Induction 
            for (int t = T - 1; t > 1; t--)
            {
                double scaler = 0;
                for (int i = 0; i < nStates; i++)
                {
                    for (int j = 0; j < nStates; j++)
                        _beta[t][i] = _transition[i][j] * _emission.EstimateLikelihood(x[t + 1], j) * _beta[t + 1][j];
                    scaler += _beta[i][t];
                }
                for (int i = 0; i < nStates; i++)
                    _beta[t][i] = _beta[t][i] / scaler;
            }
        }


        public void MaximisationStep(List<List<double>> x)
        {
            UpdateTransition();
            UpdateEmission(x);
        }

        public double EstimationStep(List<List<double>> x)
        {
            double likelihood = Forward(x);
            Backward(x);

            for (int t = 0; t < T; t++)
            {
                for (int i = 0; i < nStates; i++)
                {
                    var denominator = _alpha[t].Zip(_beta[t], (a, b) => a * b).Sum();
                    _gamma[t][i] = _alpha[t][i] * _beta[t][i] / denominator;
                    for (int j = 0; j < nStates; j++)
                        _epsilon[t][i][j] = _alpha[t][i] * _transition[i][j] * _emission.EstimateLikelihood(x[t + 1], j) *
                            _beta[t + 1][j] / denominator;
                }
            }
            return likelihood;
        }

        public void FindMaximalLikelyhood(List<List<double>> x)
        {
            double likelihoodDifferenceThreshold = 0;
            const int maxIterations = 20;
            List<double> likelihoods = new List<double>();
            double oldLikelihood = EstimationStep(x);
            likelihoods.Add(oldLikelihood);
            MaximisationStep(x);
            int iterations = 1;
            double likelihoodDifference = Double.MaxValue;
            while (iterations < maxIterations && likelihoodDifference < likelihoodDifferenceThreshold)
            {
                var newLikelihood = EstimationStep(x);
                MaximisationStep(x);
                likelihoodDifference = Math.Abs(oldLikelihood - newLikelihood);
                likelihoods.Add(newLikelihood);
                oldLikelihood = newLikelihood;
                iterations++;
            }
        }

        public List<int> BestPathViterbi(List<List<double>> x)
        {
            // Initialization 
            double[][] bestScore = CanvasCommon.Utilities.MatrixCreate(nStates, nStates);
            List<int> bestStateSequence = new List<int>(T);
            bestStateSequence.Add(0);

            // Induction 
            for (int t = 1; t < T; t++)
            {
                for (int j = 0; j < nStates; j++)
                {
                    double max = Double.MaxValue;
                    int state = 0;
                    for (int i = 0; i < nStates; i++)
                    {
                        var tmpMax = bestScore[t][i]*_transition[i][j]*_emission.EstimateLikelihood(x[t + 1], j);
                        if (tmpMax > max)
                        {
                            state = j;
                            max = tmpMax;
                        }
                    }
                    bestScore[t + 1][j] = max;
                    bestStateSequence.Add(state);
                }
            }
            return bestStateSequence;
        }
    }
}
