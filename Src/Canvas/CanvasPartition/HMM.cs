using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.Distributions;
using System.Threading.Tasks;

namespace CanvasPartition
{

    public class HiddenMarkovModel
    {
        // hidden variable parameters
        private readonly MixtureDistibution _emission;
        private readonly List<List<double>> _variances;
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


        public HiddenMarkovModel(List<List<double>> data, List<MultivariateNegativeBinomial> mixturesDistribution, List<double> haploidMeans)
        {
            // HMM set-up
            nStates = mixturesDistribution.Count;
            length = data.Count;
            _stateProbabilities = new double[nStates];
            _emission = new NegativeBinomialMixture(mixturesDistribution, haploidMeans);
            _variances = mixturesDistribution.Select(x => x.Variances).ToList();
            _transition = CanvasCommon.Utilities.MatrixCreate(nStates, nStates);
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
            _beta = CanvasCommon.Utilities.MatrixCreate(length, nStates);
            _alpha = CanvasCommon.Utilities.MatrixCreate(length, nStates);
            // marginal posterior distibution
            _gamma = CanvasCommon.Utilities.MatrixCreate(nStates, length);
            // joint posterior distibution
            _epsilon = CanvasCommon.Utilities.MatrixCreate(length, nStates, nStates);
        }

        public void AdjustTransition(int state1 = 1, int state2 = 0)
        {
            for (int i = 0; i < nStates; i++)
                _transition[state1][i] = _transition[state2][i];
        }

        public void WriteEmission()
        {
            Console.WriteLine("Emissions");
            _emission.WriteMeans();
            Console.WriteLine("Transitions");
            for (int state = 0; state < nStates; state++)
                Console.WriteLine($"Transition for state {state} = {_transition[state][state]}");
        }

        public void UpdateTransition()
        {
            double[][] numerator = CanvasCommon.Utilities.MatrixCreate(nStates, nStates);
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
            _emission.UpdateMeans(_gamma, x, _variances);
            WriteEmission();
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
                if (scaler > 0)
                {
                    for (int j = 0; j < nStates; j++)
                    {
                        _alpha[t + 1][j] = _alpha[t + 1][j] / scaler;
                    }
                }
                else
                    _alpha[t + 1][2] = 1;
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

        private void adjustTransition(double diploidTransitionProb)
        {
            if (_transition[2][2] < diploidTransitionProb)
            {
                _transition[2][2] = diploidTransitionProb;
                for (int state = 0; state < nStates; state++)
                {
                    if (state == 2)
                        continue;
                    _transition[2][state] = (1.0 - diploidTransitionProb) / (nStates - 1);
                }

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

        public void FindMaximalLikelihood(List<List<double>> x)
        {
            double likelihoodDifferenceThreshold = 0.01;
            const int maxIterations = 1;
            List<double> likelihoods = new List<double>();
            double oldLikelihood = EstimationStep(x);
            likelihoods.Add(oldLikelihood);
            MaximisationStep(x);
            int iteration = 1;
            double likelihoodDifference = Double.MaxValue;
            while (iteration < maxIterations && likelihoodDifference > likelihoodDifferenceThreshold)
            {
                var newLikelihood = EstimationStep(x);
                MaximisationStep(x);
                likelihoodDifference = Math.Abs(oldLikelihood - newLikelihood);
                likelihoods.Add(newLikelihood);
                oldLikelihood = newLikelihood;
                iteration++;
            }

            var genomeSize = Math.Pow(30, 9);
            var nonDiploidBases = 2000 * 10000;
            var diploidTransitionProb = (genomeSize - nonDiploidBases) / genomeSize;
            adjustTransition(diploidTransitionProb);
            WriteEmission();
        }

        /// <summary>
        /// BestHsmmPathViterbi helper method to calculate sojourn survival function
        /// </summary>
        public double[][] CalculateSojourn(int maxStateLength, List<int> means)
        {
            double[][] D = CanvasCommon.Utilities.MatrixCreate(maxStateLength + 1, maxStateLength + 2);
            // Store D
            for (int j = 0; j < nStates; j++)
            {
                for (int u = 1; u < maxStateLength + 1; u++)
                {
                    double x = 0;
                    for (int v = u; v < maxStateLength + 2; v++)
                        x += Math.Log(Poisson.PMF(means[j], v));
                    D[j][u] = x;
                }
                D[j][maxStateLength] = 0;
            }
            return D;
        }

        private static List<List<double>> GetStateDurationProbability(List<int> sojournMeans, int maxStateLength)
        {
            var stateDurationProbability = new List<List<double>>(sojournMeans.Count);
            var stateCounter = 0;
            foreach (var sojournMean in sojournMeans)
            {
                stateDurationProbability.Add(new List<double>(maxStateLength));
                for (var stateDuration = 0; stateDuration < maxStateLength; stateDuration++)
                    stateDurationProbability[stateCounter].Add(Math.Log(Poisson.PMF(sojournMean, stateDuration)));
                stateCounter++;
            }
            return stateDurationProbability;
        }


        /// <summary>
        /// HSMM Viterbi implementtion based on:
        /// Guedon, Y. (2003), Estimating hidden semi-Markov chains from discrete sequences, Journal of
        /// Computational and Graphical Statistics, Volume 12, Number 3, page 604-639 - 2003
        /// </summary>
        public List<int> BestHsmmPathViterbi(List<List<double>> x, List<double> haploidMeans)
        {
            // Initialization 
            var length = x.Count;
            var alpha = CanvasCommon.Utilities.MatrixCreate(nStates, length + 1);
            var bestStateDuration = new int[nStates][];
            var bestStateIndex = new int[nStates][];
            for (int i = 0; i < nStates; ++i)
            {
                bestStateIndex[i] = new int[length];
                bestStateDuration[i] = new int[length];
            }
            for (int j = 0; j < nStates; j++)
            {
                alpha[j][0] = this._stateProbabilities[j];
            }

            var maxStateLength = 110;
            var sojournMeans = new List<int>{10,10,100,50,50};
            var stateDurationProbability = GetStateDurationProbability(sojournMeans, maxStateLength);
            var sojournLastState = CalculateSojourn(maxStateLength, sojournMeans);

            double emissionSequence = 0;
            double tempEmissionSequence = 0;
            var bestState = 0;
            var firstState = true;
            var firstI = true;
            var transition = Enumerable.Repeat(1.0, nStates).ToArray();


            // Induction 
            for (int t = 1; t < length - 1; t++)
            {
                for (int j = 0; j < nStates; j++)
                {
                    emissionSequence = 0;
                    firstState = true;

                    for (int stateDuration = 1; stateDuration < Math.Min(maxStateLength, t); stateDuration++)
                    {
                        firstI = true;
                        for (int i = 0; i < nStates; i++)
                        {
                            if (i == j) continue;
                            if (Math.Log(_transition[i][j]) + alpha[i][t - stateDuration] > tempEmissionSequence || firstI)
                            {
                                tempEmissionSequence = Math.Log(_transition[i][j]) + alpha[i][t - stateDuration];
                                bestState = i;
                                firstI = false;
                            }
                        }
                        if (firstState || emissionSequence + stateDurationProbability[j][stateDuration] + tempEmissionSequence > alpha[j][t])
                        {
                            alpha[j][t] = emissionSequence + stateDurationProbability[j][stateDuration] + tempEmissionSequence;
                            bestStateDuration[j][t] = stateDuration;
                            bestStateIndex[j][t] = bestState;
                            firstState = false;
                        }
                        emissionSequence += _emission.EstimateViterbiLikelihood(x[t - stateDuration], j, haploidMeans, transition);
                    }

                    if (t + 1 <= maxStateLength)
                    {
                        if (firstState || emissionSequence + Math.Log(Poisson.PMF(sojournMeans[j], t + 1) * _stateProbabilities[j]) > alpha[j][t])
                        {
                            alpha[j][t] = emissionSequence + Math.Log(Poisson.PMF(sojournMeans[j], t + 1) * _stateProbabilities[j]);
                            bestStateDuration[j][t] = -1;
                            bestStateIndex[j][t] = -1;
                        }
                    }
                    alpha[j][t] += _emission.EstimateViterbiLikelihood(x[t], j, haploidMeans, transition);
                }
            }


            for (int j = 0; j < nStates; j++)
            {
                emissionSequence = 0;
                firstState = true;
                for (int stateDuration = 1; stateDuration < maxStateLength - 1; stateDuration++)
                {
                    firstI = true;
                    for (int i = 0; i < nStates; i++)
                    {
                        if (i == j) continue;
                        if (Math.Log(_transition[i][j]) + alpha[i][length - 1 - stateDuration] > tempEmissionSequence || firstI)
                        {
                            tempEmissionSequence = Math.Log(_transition[i][j]) + alpha[i][length - 1 - stateDuration];
                            bestState = i;
                            firstI = false;
                        }
                    }

                    if (emissionSequence + Math.Log(sojournLastState[j][Math.Min(stateDuration, maxStateLength)]) + tempEmissionSequence > alpha[j][length - 1] || firstState)
                    {
                        alpha[j][length - 1] = emissionSequence + Math.Log(sojournLastState[j][Math.Min(stateDuration, maxStateLength)]) + tempEmissionSequence;
                        bestStateDuration[j][length - 1] = stateDuration;
                        bestStateIndex[j][length - 1] = bestState;
                        firstState = false;
                    }
                    emissionSequence += _emission.EstimateViterbiLikelihood(x[length - 1 - stateDuration], j, haploidMeans, transition);
                }

                if (emissionSequence + Math.Log(sojournLastState[j][Math.Min(length - 1, maxStateLength)] * _stateProbabilities[j]) > alpha[j][length - 1] || firstState)
                {
                    alpha[j][length - 1] = emissionSequence + Math.Log(sojournLastState[j][Math.Min(length, maxStateLength)] * _stateProbabilities[j]);
                    bestStateDuration[j][length - 1] = -1;
                    bestStateIndex[j][length - 1] = -1;
                }
                alpha[j][length - 1] += _emission.EstimateViterbiLikelihood(x[length - 1], j, haploidMeans, transition);
            }

            // backtracking 
            List<int> finalStates = Enumerable.Repeat(2, length).ToList();

            int T = length - 1;
            while (bestStateIndex[bestState][T] >= 0)
            {
                for (int i = T; i >= T - bestStateDuration[bestState][T] + 1; i--)
                {
                    finalStates[i] = bestState;
                }
                var alternativeBestState = bestState;
                bestState = bestStateIndex[bestState][T];

                T -= bestStateDuration[alternativeBestState][T];
            }
            finalStates.Reverse();

            return finalStates;
        }


        public List<int> BestPathViterbi(List<List<double>> x, uint[] start, List<double> haploidMeans)
        {
            // Initialization 
            var length = x.Count;
            double[][] bestScore = CanvasCommon.Utilities.MatrixCreate(length + 1, nStates);
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
                        var tmpMax = bestScore[t - 1][i] + _emission.EstimateViterbiLikelihood(x[t], j, haploidMeans, _transition[i]);
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

            // backtracking 
            bestStates.Add(backtrack);
            while (backtrack >= 0)
            {
                bestStates.Add(bestStateSequence[backtrack][bestState]);
                backtrack--;
            }

            bestStates.Reverse();
            OutlierMask(bestStates);
            SmallSegmentsMask(bestStates);
            OversegmentationMask(bestStates);

            return bestStates;
        }

        public void OutlierMask(List<int> bestStates)
        {
            for (int k = 3; k < bestStates.Count - 3; k++)
            {
                if (bestStates[k - 2] == bestStates[k - 1] && bestStates[k + 1] == bestStates[k + 2] &&
                    bestStates[k - 1] == bestStates[k + 1] && bestStates[k] != bestStates[k - 1])
                    bestStates[k] = bestStates[k - 1];
            }

        }

        public void SmallSegmentsMask(List<int> bestStates)
        {
            List<int> bestStateDuration = Enumerable.Repeat(1, bestStates.Count + 1).ToList();
            var lastBestState = bestStates.Take(1).Single();
            int currentStateDuration = 1;
            int counter = 1;
            foreach (int bestState in bestStates.Skip(1))
            {
                if (lastBestState == bestState)
                {
                    currentStateDuration++;
                    counter++;
                    bestStateDuration[counter] = currentStateDuration;
                }
                else
                {
                    currentStateDuration = 1;
                    counter++;
                    bestStateDuration[counter] = currentStateDuration;
                    lastBestState = bestState;
                }
            }
        }

        public void OversegmentationMask(List<int> bestStates)
        {
            List<int> bestStateDuration = Enumerable.Repeat(1, bestStates.Count + 1).ToList();
            const int stateDurationCutoff = 5;
            const int halfWindow = 3;


            for (var k = halfWindow; k < bestStates.Count - halfWindow; k++)
            {
                // small segments
                if (bestStates[k] != bestStates[k + 1] && bestStates[k + 1] == bestStates[k + 2] &&
                                    bestStates[k] == bestStates[k + 2] && bestStateDuration[k] > stateDurationCutoff)
                {
                    bestStates[k + 1] = bestStates[k];
                    bestStates[k + 2] = bestStates[k];
                    bestStateDuration[k + 1] = bestStateDuration[k] + bestStateDuration[k + 1];
                    bestStateDuration[k + 2] = bestStateDuration[k + 1] + bestStateDuration[k + 2];
                }
                // oversegmentation
                if (bestStates[k] != bestStates[k + 1] && bestStates[k - 1] != bestStates[k] &&
                    bestStateDuration[k] > stateDurationCutoff)
                {
                    if (bestStates[k - 1] != bestStates[k + 1])
                    {
                        bestStates[k - 1] = bestStates[k];
                        bestStateDuration[k] = bestStateDuration[k] + bestStateDuration[k - 1];
                    }
                    else
                    {
                        bestStates[k] = bestStates[k - 1];
                        bestStateDuration[k] = bestStateDuration[k] + bestStateDuration[k - 1];
                    }
                }
            }
        }
    }
}
