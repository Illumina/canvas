using System;
using System.Collections.Generic;
using System.Text;
using CanvasCommon;

namespace CanvasPedigreeCaller
{
    public class CopyNumberModelFactory : ICopyNumberModelFactory
    {
        public ICopyNumberModel CreateModel(int numCnStates, int maxCoverage, double meanCoverage, double diploidAlleleMeanCounts)
        {

            var cnDistribution = new List<List<double>>();
            double haploidMafMean = diploidAlleleMeanCounts / 2.0;
            double haploidMean = meanCoverage / 2.0;
            double mafVariance = diploidAlleleMeanCounts * 2.5;
            double variance = meanCoverage * 2.5;
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
                cnDistribution.Add(DistributionUtilities.NegativeBinomialWrapper(haploidMean * multiplier, variance,
                    maxCoverage, adjustClumpingParameter: true));
            }

            var alleleDistribution = new Tuple<List<double>, List<double>>[numCnStates][];

            for (int i = 0; i < numCnStates; i++)
                alleleDistribution[i] = new Tuple<List<double>, List<double>>[numCnStates];

            for (int gt1 = 0; gt1 < numCnStates; gt1++)
            {
                for (int gt2 = 0; gt2 < numCnStates; gt2++)
                {
                    var gt1Probabilities = DistributionUtilities.NegativeBinomialWrapper(haploidMafMean * gt1, mafVariance, maxCoverage);
                    var gt2Probabilities = DistributionUtilities.NegativeBinomialWrapper(haploidMafMean * gt2, mafVariance, maxCoverage);
                    alleleDistribution[gt1][gt2] = new Tuple<List<double>, List<double>>(gt1Probabilities, gt2Probabilities);
                }
            }

            return new CopyNumberModel(cnDistribution, alleleDistribution, maxCoverage);
        }

    }


    public class HaplotypeCopyNumberModelFactory : ICopyNumberModelFactory
    {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="numCnStates">number of CN states considered by the model</param>
        /// <param name="maxCoverage">mean depth coverage</param>
        /// <param name="meanCoverage">max depth coverage</param>
        /// <param name="diploidAlleleMeanCounts">mean diploid allele counts</param>
        /// <returns></returns>
        public ICopyNumberModel CreateModel(int numCnStates, int maxCoverage, double meanCoverage, double diploidAlleleMeanCounts)
        {
            var cnDistribution = new List<List<double>>();
            double haploidAlleleMeanCounts = diploidAlleleMeanCounts / 2.0;
            double haploidMean = meanCoverage / 2.0;
            double mafVariance = diploidAlleleMeanCounts * 2.5;
            double variance = meanCoverage * 2.5;
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
                cnDistribution.Add(DistributionUtilities.NegativeBinomialWrapper(haploidMean * multiplier, variance,
                    maxCoverage, adjustClumpingParameter: true));
            }

            var alleleDistribution = new Tuple<List<double>, List<double>>[numCnStates][];

            for (int i = 0; i < numCnStates; i++)
                alleleDistribution[i] = new Tuple<List<double>, List<double>>[numCnStates];

            for (int gt1 = 0; gt1 < numCnStates; gt1++)
            {
                for (int gt2 = 0; gt2 < numCnStates; gt2++)
                {
                    var gt1Probabilities =
                        DistributionUtilities.NegativeBinomialWrapper(haploidAlleleMeanCounts * Math.Max(gt1, alleleStateZeroCorrector),
                            mafVariance, maxCoverage);
                    var gt2Probabilities =
                        DistributionUtilities.NegativeBinomialWrapper(haploidAlleleMeanCounts * Math.Max(gt2, alleleStateZeroCorrector),
                            mafVariance, maxCoverage);
                    alleleDistribution[gt1][gt2] =
                        new Tuple<List<double>, List<double>>(gt1Probabilities, gt2Probabilities);
                }
            }
            // meanMafCoverage * 3 will cap the coverage at 6 CN, which corresponds to 0-5 CN range captured by the model
            // this will prevent a read stacking scenario with high depth (i.e. > 1000) returning likelihoods of 0 for all models 
            int coverageCeiling = Convert.ToInt32(diploidAlleleMeanCounts * 3);
            
            // also compute distributions for total reads as a function of total copy number,
            // so we can compute likelihoods for homozygous loci
            int maxTtlCoverage = 2 * maxCoverage;
            int numTtlCnStates = 2 * numCnStates;
            var readDepthDistribution = new List<double>[numTtlCnStates];
            for (int gt1 = 0; gt1 < numTtlCnStates; gt1++)
            {
                readDepthDistribution[gt1] = DistributionUtilities.NegativeBinomialWrapper(haploidAlleleMeanCounts * gt1,
                    mafVariance, maxTtlCoverage);
            }

            return new HaplotypeCopyNumberModel(cnDistribution, alleleDistribution, coverageCeiling, readDepthDistribution, maxTtlCoverage);

        }

    }

}
