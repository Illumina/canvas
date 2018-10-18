using System;
using System.Collections.Generic;

namespace CanvasPartition
{

    public class HiddenMarkovModel
    {
        // hidden variable parameters
        private readonly MixtureDistibution _emission;
        private readonly double[][] _transition;
        private readonly double[] _stateProbabilities;
        // dimention variables
        public int nStates;
        public int length;
        public const double selfTransition = 0.99;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="data"></param>
        /// <param name="mixturesDistribution"></param>
        /// <param name="haploidMeans"></param>
        public HiddenMarkovModel(List<List<double>> data, List<MultivariateNegativeBinomial> mixturesDistribution, List<double> haploidMeans, bool isPerSample)
        {
            // HMM set-up
            nStates = mixturesDistribution.Count;
            length = data.Count;
            _stateProbabilities = new double[nStates];
            _emission = new NegativeBinomialMixture(mixturesDistribution, haploidMeans, isPerSample);
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
                _stateProbabilities[i] = 1f / nStates;
            }

            // alpha-beta algorithm
            CanvasCommon.Utilities.MatrixCreate(length, nStates);
            CanvasCommon.Utilities.MatrixCreate(length, nStates);
            // marginal posterior distibution
            CanvasCommon.Utilities.MatrixCreate(nStates, length);
            // joint posterior distibution
            CanvasCommon.Utilities.MatrixCreate(length, nStates, nStates);
        }

        /// <summary>
        /// Standard Viterbi algorithm for finding the best path through the sequence 
        /// see Rabiner, Lawrence R. "A tutorial on hidden Markov models and selected applications in speech recognition." 
        /// Proceedings of the IEEE 77.2 (1989): 257-286.
        /// </summary>
        /// <param name="depthList"></param>
        /// <param name="start"></param>
        /// <param name="haploidMeans"></param>
        /// <returns></returns>
        public List<int> BestPathViterbi(List<List<double>> depthList, uint[] start, List<double> haploidMeans)
        {
            var x = depthList;

            // Initialization 
            var size = x.Count;
            double[][] bestScore = CanvasCommon.Utilities.MatrixCreate(size, nStates);
            int[][] bestStateSequence = new int[size][];
            for (int i = 0; i < size; ++i)
                bestStateSequence[i] = new int[nStates];

            for (int j = 0; j < nStates; j++)
            {
                // This should be the score for emitting the first data element, combining the initial state prob with the emission prob.
                // The right way to make the change here is to refactor such that we can get the emission probability separate from the transition probability,
                // but for now just hack: subtract off the transition probability
                bestScore[0][j] = Math.Log(this._stateProbabilities[j]) + _emission.EstimateViterbiLikelihood(x[0], j, haploidMeans, _transition[0]) - Math.Log(_transition[0][j]);
                bestStateSequence[0][j] = -1;
            }

            // Induction 
            for (int t = 1; t < size; t++)
            {
                for (int j = 0; j < nStates; j++)
                {
                    int state = 0;
                    double max = Double.MinValue;
                    for (int i = 0; i < nStates; i++)
                    {
                        var vitLogL = _emission.EstimateViterbiLikelihood(x[t], j, haploidMeans, _transition[i]);
                        var tmpMax = bestScore[t - 1][i] + vitLogL;
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

            int bestState = -1;
            var bestStates = new List<int>(size);
            var max1 = Double.MinValue;
            for (int i = 0; i < nStates; i++)
            {

                var tmpMax = bestScore[size - 1][i];
                if (tmpMax > max1)
                {
                    bestState = i;
                    max1 = tmpMax;
                }
            }

            // backtracking 
            var backtrack = size - 1;
            while (backtrack > 0)
            {
                bestStates.Add(bestState);
                bestState = bestStateSequence[backtrack][bestState];
                backtrack--;
            }
            bestStates.Add(bestState);

            bestStates.Reverse();
            return bestStates;
        }
    }
}
