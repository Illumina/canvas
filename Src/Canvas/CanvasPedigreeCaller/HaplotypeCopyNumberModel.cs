using System;
using System.Collections.Generic;
using System.Linq;
using CanvasCommon;

namespace CanvasPedigreeCaller
{

    public class HaplotypeCopyNumberModel : ICopyNumberModel
    {
        private readonly List<List<double>> _cnDistribution;
        private readonly Tuple<List<double>, List<double>>[][] _alleleDistribution;
        private readonly int _maxCoverage;

        public HaplotypeCopyNumberModel(List<List<double>> cnDistribution, Tuple<List<double>, List<double>>[][] alleleDistribution, int maxCoverage)
        {
            _cnDistribution = cnDistribution;
            _alleleDistribution = alleleDistribution;
            _maxCoverage = maxCoverage;
        }

        public double GetTotalCopyNumberLikelihoods(double segmentMedianBinCoverage, Genotype totalCopyNumberGenotype)
        {
            return _cnDistribution.Select(x => x[Convert.ToInt32(segmentMedianBinCoverage)]).ElementAt(totalCopyNumberGenotype.TotalCopyNumber);
        }

        public double GetGenotypeLikelihood(Balleles gtObservedCounts, PhasedGenotype gtModelCount)
        {
            {
                const double lohRefModelPenaltyTerm = 0.5;
                double currentLikelihood = 0;
                foreach (var gtCount in gtObservedCounts.GetAlleleCounts())
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

