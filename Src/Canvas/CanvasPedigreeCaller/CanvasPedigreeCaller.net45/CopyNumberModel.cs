using System;
using System.Collections.Generic;
using System.Linq;
using CanvasCommon;

namespace CanvasPedigreeCaller
{
    public class CopyNumberModel
    {
        public List<Tuple<int,int>> Genotypes = new List<Tuple<int, int>>();
        public List<List<double>> CnDistribution = new List<List<double>>();
        readonly Tuple<List<double>, List<double>> [][] _alleleDistribution;

        public CopyNumberModel(int numCnStates, double haploidMean, double haploidMafMean, double variance, double mafVariance, int maxValue)
        {
            for (int copyNumber = 0; copyNumber  < numCnStates; copyNumber ++)
            {
                var multiplier = 1.0;
                if (copyNumber == 1)
                    multiplier = 0.9;
                if (copyNumber == 3)
                    multiplier = 1.1;
                CnDistribution.Add(DistributionUtilities.NegativeBinomialWrapper(haploidMean * copyNumber * multiplier, variance, maxValue, adjustR: true));
            }

            _alleleDistribution = new Tuple<List<double>, List<double>>[numCnStates][];

            for (int i = 0; i < numCnStates; i++)
                _alleleDistribution[i] = new Tuple<List<double>, List<double>>[numCnStates];

            for (int gt1 = 0; gt1 < numCnStates; gt1++)
            {
                for (int gt2 = 0; gt2 < numCnStates; gt2++)
                {
                    var gt1Probabilities = DistributionUtilities.NegativeBinomialWrapper(haploidMafMean * gt1, mafVariance, maxValue);
                    var gt2Probabilities = DistributionUtilities.NegativeBinomialWrapper(haploidMafMean * gt2, mafVariance, maxValue);
                    _alleleDistribution[gt1][gt2] = new Tuple<List<double>, List<double>>(gt1Probabilities, gt2Probabilities);
                }
            }
        }


        public List<double> GetCnLikelihood(double dimension)
        {
            return CnDistribution.Select(x => x[Convert.ToInt32(dimension)]).ToList();
        }

        public double[][] GetMedianGtLikelihood(Tuple<int,int> gtCounts)
        {
            int nrows = _alleleDistribution.Length;
            int ncols = _alleleDistribution.First().Length;
            double[][] likelihood = Utilities.MatrixCreate(nrows, ncols);

            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    if (_alleleDistribution[i][j] != null)
                        likelihood[i][j] = _alleleDistribution[i][j].Item1[gtCounts.Item1]*
                                           _alleleDistribution[i][j].Item2[gtCounts.Item2];
                    else
                        likelihood[i][j] = 0;
                }
            }
            return likelihood;
        }

        public double GetGtLikelihood(List<Tuple<int, int>> gtCounts, List<Tuple<int, int>> gtModelCounts, ref int selectedGtState, int maxCoverage)
        {
            double maximalLikelihood = 0;
            var bestGtStates = Enumerable.Repeat(0, gtModelCounts.Count).ToList();
            foreach (var gtCount in gtCounts)
            {
                double currentLikelihood = 0;
                var bestLikelihood = currentLikelihood;
                var bestGtState = 0;
                var gtCounter = 0;
                foreach (var gtModelCount in gtModelCounts)
                {
                    var rowId = Math.Min(gtCount.Item1, maxCoverage - 1);
                    var colId = Math.Min(gtCount.Item2, maxCoverage - 1);
                    currentLikelihood = _alleleDistribution[gtModelCount.Item1][gtModelCount.Item2].Item1[rowId] *
                     _alleleDistribution[gtModelCount.Item1][gtModelCount.Item2].Item2[colId];

                    if (currentLikelihood > bestLikelihood)
                    {
                        bestLikelihood = currentLikelihood;
                        bestGtState = gtCounter;
                    }
                    gtCounter++;
                }
                maximalLikelihood += bestLikelihood;
                bestGtStates[bestGtState] += 1;
            }
            selectedGtState = bestGtStates.IndexOf(bestGtStates.Max());
            return maximalLikelihood;
        }
    }
}
