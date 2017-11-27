using System;
using System.Collections.Generic;
using System.Linq;
using CanvasCommon;

namespace CanvasPedigreeCaller
{

    public class HaplotypeCopyNumberModel : ICopyNumberModel
    {
        private readonly List<List<double>> _cnDistribution = new List<List<double>>();
        private Tuple<List<double>, List<double>>[][] _alleleDistribution;
        private int _maxCoverage;


        public void InitializeModel(int numCnStates, int maxCoverage, double meanCoverage, double meanMafCoverage)
        {
            double haploidMafMean = meanMafCoverage / 2.0;
            double haploidMean = meanCoverage / 2.0;
            double mafVariance = meanMafCoverage * 2.5;
            double variance = meanCoverage * 2.5;
            _maxCoverage = maxCoverage;
            const double alleleStateZeroCorrector = 0.1; // prevent from using zero as a mean of model distribution
            for (int copyNumber = 0; copyNumber < numCnStates; copyNumber++)
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
                _cnDistribution.Add(DistributionUtilities.NegativeBinomialWrapper(haploidMean * multiplier, variance,
                    maxCoverage,
                    adjustClumpingParameter: true));
            }

            _alleleDistribution = new Tuple<List<double>, List<double>>[numCnStates][];

            for (int i = 0; i < numCnStates; i++)
                _alleleDistribution[i] = new Tuple<List<double>, List<double>>[numCnStates];

            for (int gt1 = 0; gt1 < numCnStates; gt1++)
            {
                for (int gt2 = 0; gt2 < numCnStates; gt2++)
                {
                    var gt1Probabilities =
                        DistributionUtilities.NegativeBinomialWrapper(haploidMafMean * Math.Max(gt1, alleleStateZeroCorrector),
                            mafVariance, maxCoverage);
                    var gt2Probabilities =
                        DistributionUtilities.NegativeBinomialWrapper(haploidMafMean * Math.Max(gt2, alleleStateZeroCorrector),
                            mafVariance, maxCoverage);
                    _alleleDistribution[gt1][gt2] =
                        new Tuple<List<double>, List<double>>(gt1Probabilities, gt2Probabilities);
                }
            }
        }

        public double GetGtLikelihoodScore(List<Tuple<int, int>> gtObservedCounts, List<PhasedGenotype> gtModelCounts, ref int? selectedGtState)
        {
            throw new NotImplementedException();
        }

        public List<double> GetCnLikelihood(double dimension)
        {
            return _cnDistribution.Select(x => x[Convert.ToInt32(dimension)]).ToList();
        }

        public double GetGtLikelihood(List<Tuple<int, int>> gtObservedCounts, PhasedGenotype gtModelCount)
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

