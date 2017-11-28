using System;
using System.Collections.Generic;
using System.Text;
using CanvasCommon;

namespace CanvasPedigreeCaller
{
    public class CopyNumberModelFactory : ICopyNumberModelFactory
    {
        public ICopyNumberModel CreateModel(int numCnStates, int maxCoverage, double meanCoverage, double meanMafCoverage)
        {

            var cnDistribution = new List<List<double>>();
            double haploidMafMean = meanMafCoverage / 2.0;
            double haploidMean = meanCoverage / 2.0;
            double mafVariance = meanMafCoverage * 2.5;
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
                    var gt2Probabilities = DistributionUtilities.NegativeBinomialWrapper(haploidMafMean * gt1, mafVariance, maxCoverage);
                    alleleDistribution[gt1][gt2] = new Tuple<List<double>, List<double>>(gt1Probabilities, gt2Probabilities);
                }
            }

            return new CopyNumberModel(cnDistribution, alleleDistribution, maxCoverage);
        }

    }

    public class HaplotypeCopyNumberModelFactory : ICopyNumberModelFactory
    {
        public ICopyNumberModel CreateModel(int numCnStates, int maxCoverage, double meanCoverage, double meanMafCoverage)
        {
            var cnDistribution = new List<List<double>>();
            double haploidMafMean = meanMafCoverage / 2.0;
            double haploidMean = meanCoverage / 2.0;
            double mafVariance = meanMafCoverage * 2.5;
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
                    maxCoverage,
                    adjustClumpingParameter: true));
            }

            var alleleDistribution = new Tuple<List<double>, List<double>>[numCnStates][];

            for (int i = 0; i < numCnStates; i++)
                alleleDistribution[i] = new Tuple<List<double>, List<double>>[numCnStates];

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
                    alleleDistribution[gt1][gt2] =
                        new Tuple<List<double>, List<double>>(gt1Probabilities, gt2Probabilities);
                }
            }

            return new HaplotypeCopyNumberModel(cnDistribution, alleleDistribution, maxCoverage);
        }

    }

}
