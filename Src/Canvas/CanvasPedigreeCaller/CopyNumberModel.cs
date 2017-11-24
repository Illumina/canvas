using System;
using System.Collections.Generic;
using System.Linq;
using CanvasCommon;

namespace CanvasPedigreeCaller
{
    public class CopyNumberModel
    {
        private readonly List<List<double>> _cnDistribution = new List<List<double>>();
        private readonly Tuple<List<double>, List<double>> [][] _alleleDistribution;
        private readonly int _maxCoverage;

        public CopyNumberModel(int numCnStates, double meanMafCoverage, double meanCoverage, int maxCoverage)
        {
            double haploidMafMean = meanMafCoverage / 2.0;
            double haploidMean = meanCoverage / 2.0;
            double mafVariance = meanMafCoverage * 2.5;
            double variance = meanCoverage * 2.5;
            _maxCoverage = maxCoverage;

            for (int copyNumber = 0; copyNumber  < numCnStates; copyNumber ++)
            {
                var multiplier = copyNumber * 1.0;
                // lower haploid mean by 10% to offset FP CN=1 calls 
                if (copyNumber == 1)
                    multiplier *= 0.9;
                // allow for low bin counts for CN=0, i.e. due to imprecise breakpoints  
                if (copyNumber == 0)
                    multiplier = 0.1;
                // increase triploid mean by 10% to offset FP CN=3 calls 
                if (copyNumber == 3)
                    multiplier *= 1.1;
                _cnDistribution.Add(DistributionUtilities.NegativeBinomialWrapper(haploidMean * multiplier, variance, maxCoverage, 
                    adjustClumpingParameter: true));
            }

            _alleleDistribution = new Tuple<List<double>, List<double>>[numCnStates][];

            for (int i = 0; i < numCnStates; i++)
                _alleleDistribution[i] = new Tuple<List<double>, List<double>>[numCnStates];

            for (int gt1 = 0; gt1 < numCnStates; gt1++)
            {
                for (int gt2 = 0; gt2 < numCnStates; gt2++)
                {
                    var gt1Probabilities = DistributionUtilities.NegativeBinomialWrapper(haploidMafMean * gt1, mafVariance, maxCoverage);
                    var gt2Probabilities = DistributionUtilities.NegativeBinomialWrapper(haploidMafMean * gt2, mafVariance, maxCoverage);
                    _alleleDistribution[gt1][gt2] = new Tuple<List<double>, List<double>>(gt1Probabilities, gt2Probabilities);
                }
            }
        }


        public List<double> GetCnLikelihood(double dimension)
        {
            return _cnDistribution.Select(x => x[Convert.ToInt32(dimension)]).ToList();
        }

        public double GetGtLikelihoodScore(List<Tuple<int, int>> gtObservedCounts, List<PhasedGenotype> gtModelCounts, ref int? selectedGtState)
        {
            const int maxGQscore = 60;
            var gtLikelihoods = Enumerable.Repeat(0.0, gtModelCounts.Count).ToList();
            var gtModelCounter = 0;
            foreach (var gtModelCount in gtModelCounts)
            {
                gtLikelihoods[gtModelCounter] = GetCurrentGtLikelihood(gtObservedCounts, gtModelCount);
                gtModelCounter++;
            }
            if (!selectedGtState.HasValue)
                selectedGtState = gtLikelihoods.IndexOf(gtLikelihoods.Max());
            double normalizationConstant = gtLikelihoods.Sum();
            double gqscore = -10.0 * Math.Log10((normalizationConstant - gtLikelihoods[selectedGtState.Value]) / normalizationConstant);
            if (Double.IsInfinity(gqscore) | gqscore > maxGQscore)
                gqscore = maxGQscore;
            return Double.IsNaN(gqscore) || Double.IsInfinity(gqscore) ? 0 : gqscore;
        }

        public double GetCurrentGtLikelihood(List<Tuple<int, int>> gtObservedCounts, PhasedGenotype gtModelCount)
        {
            {
                const double lohRefModelPenaltyTerm = 0.5;
                double currentLikelihood = 0;
                foreach (var gtCount in gtObservedCounts)
                {
                    int rowId = Math.Min(gtCount.Item1, _maxCoverage - 1);
                    int colId = Math.Min(gtCount.Item2, _maxCoverage - 1);
                    if (gtModelCount.CopyNumberA == 1 && gtModelCount.CopyNumberB == 1)
                        currentLikelihood += new List<double>
                        {
                            _alleleDistribution[1][1].Item1[rowId] *
                            _alleleDistribution[1][1].Item2[colId],
                            _alleleDistribution[2][0].Item1[rowId] *
                            _alleleDistribution[2][0].Item2[colId],
                            _alleleDistribution[0][2].Item1[rowId] *
                            _alleleDistribution[0][2].Item2[colId]
                        }.Max() * lohRefModelPenaltyTerm;
                    else
                        currentLikelihood += _alleleDistribution[gtModelCount.CopyNumberA][gtModelCount.CopyNumberB]
                                                 .Item1[rowId] *
                                             _alleleDistribution[gtModelCount.CopyNumberA][gtModelCount.CopyNumberB]
                                                 .Item2[colId];
                }
                return currentLikelihood;
            }
        }
    }
}

