using System;
using System.Collections.Generic;
using System.Linq;
using CanvasCommon;

namespace CanvasPedigreeCaller
{

    public class CopyNumberModel : ICopyNumberModel
    {
        private readonly List<List<double>> _cnDistribution;
        private readonly Tuple<List<double>, List<double>>[][] _alleleDistribution;
        private readonly int _maxCoverage;

        public CopyNumberModel(List<List<double>> cnDistribution, Tuple<List<double>, List<double>>[][] alleleDistribution, int maxCoverage)
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
            double currentLikelihood = 0;
            foreach (var gtCount in gtObservedCounts.GetAlleleCounts())
            {
                int rowId = Math.Min(gtCount.Item1, _maxCoverage - 1);
                int colId = Math.Min(gtCount.Item2, _maxCoverage - 1);
                currentLikelihood += Math.Log10(_alleleDistribution[gtModelCount.CopyNumberA][gtModelCount.CopyNumberB].Item1[rowId] *
                                     _alleleDistribution[gtModelCount.CopyNumberA][gtModelCount.CopyNumberB].Item2[colId]);
            }
            return currentLikelihood;
        }
    }
}

