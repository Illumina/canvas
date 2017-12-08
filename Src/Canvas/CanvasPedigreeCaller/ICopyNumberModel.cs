using CanvasCommon;

namespace CanvasPedigreeCaller
{
    public interface ICopyNumberModel
    {
        double GetGenotypeLikelihood(Balleles segmentAlleleReadCounts, PhasedGenotype phasedGenotype);
        double GetTotalCopyNumberLikelihoods(double segmentMedianBinCoverage, Genotype totalCopyNumberGenotype);
    }

    public interface ICopyNumberModelFactory
    {
        ICopyNumberModel CreateModel(int numCnStates, int maxCoverage, double haploidMean, double haploidMafMean);
    }
  
}