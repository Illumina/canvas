using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Principal;
using CanvasCommon;
using MathNet.Numerics.Distributions;
using System.Threading;
using System.Threading.Tasks;
using Isas.Shared.Utilities;

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
        public double[][] Covariance   // the Mean property
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
            int stateCounter = -1;
            foreach (MultivariateGaussianDistribution gaussianDistribution in _gaussianDistributions)
            {
                stateCounter++;
                if (stateCounter == 1)
                    continue;
                gaussianDistribution.UpdateMean(gamma[stateCounter], x);
            }
        }
        public void WriteMeans()
        {
            foreach (MultivariateGaussianDistribution gaussianDistribution in _gaussianDistributions)
                for (int i = 0; i < gaussianDistribution.Mean.Length; i++)
                    Console.WriteLine($"Gaussian mean for state {i} = {gaussianDistribution.Mean[i]}");
        }
        public void UpdateCovariances(double[][] gamma, List<List<double>> x)
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
        public int length;
        public const double selfTransition = 0.99;


        public HiddenMarkovModel(List<List<double>> data, List<MultivariateGaussianDistribution> gaussianMixtures)
        {
            // HMM set-up
            nStates = gaussianMixtures.Count;
            length = data.Count;
            _stateProbabilities = new double[nStates];
            _emission = new GaussianMixture(gaussianMixtures);
            _transition = MatrixUtilities.MatrixCreate(nStates, nStates);
            for (int i = 0; i < nStates; i++)
            {
                for (int j = 0; j < nStates; j++)
                {
                    if (i == j)
                        _transition[i][j] = selfTransition;
                    else
                        _transition[i][j] = (1.0 - selfTransition) / (nStates - 1);
                }
            }

            // alpha-beta algorithm
            _beta = MatrixUtilities.MatrixCreate(length, nStates);
            _alpha = MatrixUtilities.MatrixCreate(length, nStates);
            // marginal posterior distibution
            _gamma = MatrixUtilities.MatrixCreate(nStates, length);
            // joint posterior distibution
            _epsilon = MatrixUtilities.MatrixCreate(length, nStates, nStates);
        }

        public void AdjustTransition(int state1 = 1, int state2 = 0)
        {
            for (int i = 0; i < nStates; i++)
                _transition[state1][i] = _transition[state2][i];
        }

        public void writeEmission()
        {
            Console.WriteLine("Emissions");
            _emission.WriteMeans();
            Console.WriteLine("Transitions");
            for (int state = 0; state < nStates; state++)
                Console.WriteLine($"Transitin for state {state} = {_transition[state][state]}");
        }

        public void UpdateTransition()
        {
            double[][] numerator = MatrixUtilities.MatrixCreate(nStates, nStates);
            double[] denominator = new double[nStates];
            for (int t = 0; t < length - 1; t++)
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
                _alpha[0][j] = 1 / (double)nStates;

            // Induction 
            for (int t = 0; t < length - 1; t++)
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
            return likelihood;
        }
        public void Backward(List<List<double>> x)
        {
            // Initialization 
            for (int j = 0; j < nStates; j++)
                _beta[length - 1][j] = 1 / (double)nStates;


            // Induction 
            int t = length - 2;
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
            double likelihood = 0;
            Parallel.Invoke(() => likelihood = Forward(x), () => Backward(x));

            for (int t = 0; t < length - 1; t++)
            {
                for (int i = 0; i < nStates; i++)
                {
                    _gamma[i][t] = 0;

                    for (int j = 0; j < nStates; j++)
                    {
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
            const int maxIterations = 3;
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
                Console.WriteLine($"{DateTime.Now} Parameters for chromosome {chr}");
                writeEmission();
            }
            Console.WriteLine($"Completed EM step for chromosome {chr}");
        }

        public double[][] CalcStoreD(int M, List<int> means)
        {
            double[][] D = CanvasCommon.MatrixUtilities.MatrixCreate(M, M+1);
            // Store D
            for (int j = 0; j < nStates; j++)
            {
                for (int u = 1; u < M; u++)
                {
                    double x = 0;
                    for (int v = u; v < M + 1; v++)
                        x += Math.Log(Poisson.PMF(means[j], v));
                    D[j][u] = x;
                }
                    D[j][M] = 0;
            }
            return D;
        }



        public List<int> BestHsmmPathViterbi(List<List<double>> x)
        {

            // Initialization 
            AdjustTransition();
            var length = x.Count;
            double[][] bestScore = MatrixUtilities.MatrixCreate(length + 1, nStates);
            int[][] maxU = new int[length + 1][];
            int[][] maxI = new int[length + 1][];
            int[][] bestStateSequence = new int[length + 1][];
            for (int i = 0; i < length; ++i)
            {
                bestStateSequence[i] = new int[nStates];
                maxI[i] = new int[nStates];
                maxU[i] = new int[nStates];

            }

            int maxStateLength = 10000;
            int cneq2Length = 5000;
            int cnnq2Length = 10;

            List<int> means = new List<int>(nStates);
            for (int i = 0; i < nStates; i++)
                means.Add(i == 2 ? cneq2Length : cnnq2Length);

            double[][] D = CalcStoreD(maxStateLength, means);
            double emissionSequence = 0;
            double tempEmissionSequence = 0;
            int bestState = 0;
            var firstState = true;
            var firstI = true;

            // Induction 
            for (int t = 0; t < length - 1; t++)
            {
                for (int j = 0; j < nStates; j++)
                {
                    emissionSequence = 0;
                    firstState= true;

                    for (int u = 1; u < Math.Min(maxStateLength, length); u++)
                    {
                        firstI = true;
                        for (int i = 0; i <= nStates; i++)
                        {
                            if (i != j)
                            {
                                if (Math.Log(_transition[i][j]) + bestScore[i][t - u] > tempEmissionSequence || firstI)
                                {
                                    tempEmissionSequence = Math.Log(_transition[i][j]) + bestScore[i][t - u];
                                    bestState = i;
                                    firstI = false;
                                }
                            }
                        }
                        if (firstState ||
                            (emissionSequence + Math.Log(Poisson.PMF(means[j], u)) + tempEmissionSequence >
                             bestScore[j][t]))
                        {
                            bestScore[j][t] = emissionSequence + Math.Log(Poisson.PMF(means[j], u)) +
                                              tempEmissionSequence;
                            maxU[j][t] = u;
                            maxI[j][t] = bestState;
                            firstState = false;
                        }
                        emissionSequence += Math.Log(_emission.EstimateLikelihood(x[t - u], j));
                    }

                    if (t + 1 <= maxStateLength)
                    {
                        if (firstState ||
                            (emissionSequence + Math.Log(Poisson.PMF(means[j], t + 1) * _stateProbabilities[j]) > bestScore[j][t]))
                        {
                            bestScore[j][t] = emissionSequence + Math.Log(Poisson.PMF(means[j], t + 1) * _stateProbabilities[j]);
                            maxU[j][t] = -1;
                            maxI[j][t] = -1;
                        }
                    }
                    bestScore[j][t] += Math.Log(_emission.EstimateLikelihood(x[t], j));
                    }
                }


            for (int j = 0; j < nStates; j++)
            {
                emissionSequence = 0;
                firstState = true;
                for (int u = 1; u <= length - 1; u++)
                {
                    firstI = true;
                    for (int i = 0; i < nStates; i++)
                    {
                        if (i != j)
                            if ((Math.Log(_transition[i][j]) + bestScore[i][length - 1 - u] > tempEmissionSequence) || firstI)
                            {
                                tempEmissionSequence = Math.Log(_transition[i][j]) + bestScore[i][length - 1 - u];
                                bestState = i;
                                firstI = false;
                            }
                    }

                    if ((emissionSequence + Math.Log(D[j][Math.Max(u,maxStateLength)]) + tempEmissionSequence > bestScore[j][length - 1]) || firstState)
                    {
                        bestScore[j][length - 1] = emissionSequence + Math.Log(D[j][Math.Max(u, maxStateLength)]) + tempEmissionSequence;
                        maxU[j][length - 1] = u;
                        maxI[j][length - 1] = bestState;
                        firstState = false;
                    }
                    emissionSequence += Math.Log(_emission.EstimateLikelihood(x[length - 1 - u], j));
                }

                if ((emissionSequence + Math.Log(D[j][Math.Max(length - 1, maxStateLength)] * _stateProbabilities[j]) > bestScore[j][length - 1]) || firstState)
                {
                    bestScore[j][length - 1] = emissionSequence + Math.Log(D[j][Math.Max(length, maxStateLength)] * _stateProbabilities[j]);
                    maxU[j][length - 1] = -1;
                    maxI[j][length - 1] = -1;
                }
                bestScore[j][length - 1] += Math.Log(_emission.EstimateLikelihood(x[length - 1], j));
            }
            // Save result
            List<int> bestStates = new List<int>(length);

            int T = length - 1;
            while (maxI[bestState][T] >= 0)
            {
                for (int i = T; i >= T - maxU[bestState][T] + 1; i--)
                {
                    bestStates[i] = bestState;
                }
                var alternativeBestState = bestState;
                bestState = maxI[bestState][T];

                T -= maxU[alternativeBestState][T];
            }
            bestStates.Reverse();

            return bestStates;
        }


        public List<int> BestPathViterbi(List<List<double>> x)
        {
            // Initialization 
            AdjustTransition();
            var length = x.Count;
            double[][] bestScore = MatrixUtilities.MatrixCreate(length + 1, nStates);
            int[][] bestStateSequence = new int[length + 1][];
            for (int i = 0; i < length; ++i)
                bestStateSequence[i] = new int[nStates];

            for (int j = 0; j < nStates; j++)
            {
                bestScore[0][j] = this._stateProbabilities[j];
                bestStateSequence[0][j] = 0;
            }

            // Induction 
            for (int t = 1; t < length - 1; t++)
            {
                for (int j = 0; j < nStates; j++)
                {
                    int state = 0;
                    double max = Double.MinValue;
                    for (int i = 0; i < nStates; i++)
                    {
                        var tmpMax = bestScore[t - 1][i] + Math.Log(_transition[i][j]) + Math.Log(_emission.EstimateLikelihood(x[t], j));
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

            var backtrack = length - 1;
            int bestState = 0;
            List<int> bestStates = new List<int>(length + 1);
            for (int i = 0; i < nStates; i++)
            {
                double max = Double.MinValue;
                var tmpMax = bestScore[length][i];
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
