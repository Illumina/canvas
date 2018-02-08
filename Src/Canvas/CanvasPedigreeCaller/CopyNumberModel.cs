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

        public double GetGenotypeLogLikelihood(Balleles gtObservedCounts, PhasedGenotype gtModelCount)
        {
            double currentLikelihood = 0;
            double minLogLikelihood = Math.Log(1.0/Double.MaxValue);
            foreach (var gtCount in gtObservedCounts.GetAlleleCounts())
            {
                int rowId = Math.Min(gtCount.Item1, _maxCoverage - 1);
                int colId = Math.Min(gtCount.Item2, _maxCoverage - 1);
                // CopyNumberModel does not account for allele CN, aggregate into MCC
                double aAllele = Math.Log(_alleleDistribution[gtModelCount.CopyNumberA][gtModelCount.CopyNumberB]
                                                  .Item1[rowId] * _alleleDistribution[gtModelCount.CopyNumberA][gtModelCount.CopyNumberB].Item2[colId]);
                double bAllele = Math.Log(_alleleDistribution[gtModelCount.CopyNumberB][gtModelCount.CopyNumberA]
                                                  .Item1[rowId] * _alleleDistribution[gtModelCount.CopyNumberB][gtModelCount.CopyNumberA].Item2[colId]);
                currentLikelihood += Math.Max(minLogLikelihood, Math.Max(aAllele, bAllele));
                
            }
            return currentLikelihood;
        }
    }
}

