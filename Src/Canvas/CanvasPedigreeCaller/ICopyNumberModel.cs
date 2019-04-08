using CanvasCommon;

namespace CanvasPedigreeCaller
{
    public interface ICopyNumberModel
    {
        double GetGenotypeLogLikelihood(Balleles segmentAlleleReadCounts, PhasedGenotype phasedGenotype);
        double GetTotalCopyNumberLikelihoods(double segmentMedianBinCoverage, Genotype totalCopyNumberGenotype);
        int GetCoverageBound();
        int GetMaxCopyNumber();
    }

    public interface ICopyNumberModelFactory
    {
        ICopyNumberModel CreateModel(int numCnStates, int maxCoverage, double haploidMean, double haploidMafMean);
    }
  
}