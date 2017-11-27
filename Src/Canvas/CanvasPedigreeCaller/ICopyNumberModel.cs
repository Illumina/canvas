using System;
using System.Collections.Generic;
using CanvasCommon;
using Isas.Framework.DataTypes.Maps;

namespace CanvasPedigreeCaller
{
    public interface ICopyNumberModel
    {
        double GetGtLikelihood(List<Tuple<int, int>> gtObservedCounts, PhasedGenotype gtModelCount);
        List<double> GetCnLikelihood(double dimension);
        void InitializeModel(int numCnStates, int maxCoverage, double haploidMean, double variance, double haploidMafMean, double mafVariance);

        double GetGtLikelihoodScore(List<Tuple<int, int>> gtObservedCounts, List<PhasedGenotype> gtModelCounts, ref int? selectedGtState);
    }
}