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
        private readonly List<double>[] _totalAlleleCountsDistribution;
        private readonly int _maxAlleleCounts;
        private readonly List<double> _logFactorial;


        public HaplotypeCopyNumberModel(List<List<double>> cnDistribution, Tuple<List<double>, List<double>>[][] alleleDistribution,
            int maxCoverage, List<double>[] totalAlleleCountsDistribution, int maxAlleleCounts)
        {
            _cnDistribution = cnDistribution;
            _alleleDistribution = alleleDistribution;
            _maxCoverage = maxCoverage;
            _totalAlleleCountsDistribution = totalAlleleCountsDistribution;
            _maxAlleleCounts = maxAlleleCounts;

            _logFactorial = new List<double>();
            _logFactorial.Add(0.0);
            _logFactorial.Add(0.0);
            for (int i = 2; i <= _maxCoverage * 2; i++)
            {
                _logFactorial.Add(Math.Log(i));
            }
        }

        public double GetTotalCopyNumberLikelihoods(double segmentMedianBinCoverage, Genotype totalCopyNumberGenotype)
        {
            return _cnDistribution.Select(x => x[Convert.ToInt32(segmentMedianBinCoverage)]).ElementAt(totalCopyNumberGenotype.TotalCopyNumber);
        }

        public double GetGenotypeLogLikelihood(Balleles gtObservedCounts, PhasedGenotype gtModelCount)
        {
                double minLogLikelihood = Math.Log(1.0 / Double.MaxValue);
                double currentLogLikelihood = 0;
                foreach (var gtCount in gtObservedCounts.GetAlleleCounts())
                {
                    int rowId = Math.Min(gtCount.Item1, _maxCoverage - 1);
                    int colId = Math.Min(gtCount.Item2, _maxCoverage - 1);
                    int numHapsNonZero = (gtModelCount.CopyNumberA > 0 ? 1 : 0) + (gtModelCount.CopyNumberB > 0 ? 1 : 0);
                    double likelihoodThisLocus = 0;
                    // the observations can arise from a het locus, if both copy numbers are positive
                    if (numHapsNonZero == 2)
                    {
                        // Given a variant locus with two haplotypes, we have a roughly 2/3 chance of it being het.
                        // Alleles have 50:50 chance of being 'A' or 'B'.
                        // We ignore error terms, as they should have a negligible impact here.
                        likelihoodThisLocus += 1.0 / 3.0 *
                                               (
                                                   _alleleDistribution[gtModelCount.CopyNumberA][gtModelCount.CopyNumberB].Item1[rowId] *
                                                   _alleleDistribution[gtModelCount.CopyNumberA][gtModelCount.CopyNumberB].Item2[colId]
                                                   +
                                                   _alleleDistribution[gtModelCount.CopyNumberA][gtModelCount.CopyNumberB].Item1[colId] *
                                                   _alleleDistribution[gtModelCount.CopyNumberA][gtModelCount.CopyNumberB].Item2[rowId]
                                               );
                    }
                    // they can also arise from a hom locus in various ways
                    if (numHapsNonZero > 0)
                    {
                        // these should be constants to avoid calling Log over and over.
                        double logErrorProb = Math.Log(0.01);
                        double logNoErrorProb = Math.Log(.99);
                        // If both haplotypes have non-zero depth and the locus is non-ref, a locus has a prior prob of 1/3 of being hom,
                        // assuming a well-mixed population.  We could adjust for observed het:hom, but we do not at this time.
                        // Of course, if only one haplotype has non-zero depth, it must be hom.
                        double priorFactorHom = numHapsNonZero == 2 ? 1.0 / 3.0 : 1.0;
                        // limit ttlReads to maxTotalDepth as that is all we have _readDepth probabilities for
                        int ttlReads = Math.Min(rowId + colId, _maxAlleleCounts);
                        int ttlCN = gtModelCount.CopyNumberA + gtModelCount.CopyNumberB;
                        // Split the likelihood into two parts:
                        // First, the probability of getting the observed total number of reads, given the total copy number
                        double probTtlReadDepth = _totalAlleleCountsDistribution[ttlCN][ttlReads];
                        // Second, the probability of the observed per-allele read counts assuming one of the alleles is an error.
                        // The calculation here is simply binomial, in log space
                        double logProbCountAErrors = LogCombinations(rowId, colId) + rowId * logErrorProb + colId * logNoErrorProb;
                        double logProbCountBErrors = LogCombinations(rowId, colId) + colId * logErrorProb + rowId * logNoErrorProb;

                        likelihoodThisLocus += priorFactorHom * probTtlReadDepth * (
                                                   Math.Exp(logProbCountAErrors) + Math.Exp(logProbCountBErrors));
                    }
                    else
                    {
                        // if copy number is 0, any reads must be mismappings, all bets are off ... just return a constant;
                        // that constant might as well be 1 -- returning 0
                        likelihoodThisLocus = 1;
                    }

                    likelihoodThisLocus = Math.Max(minLogLikelihood, likelihoodThisLocus);
                    currentLogLikelihood += Math.Log(likelihoodThisLocus);
            }
            return currentLogLikelihood;
        }
        private double LogCombinations(int a, int b)
        {
            // N.B. This is defined as the combinations of a red plus b black, rather than 
            // the more-typical definitions as k red of n (red+black).
            // Also note: LogCombinations(a,b) == LogCombinations(b,a)
            return _logFactorial[a + b] - _logFactorial[a] - _logFactorial[b];
        }
    }
}

