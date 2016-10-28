using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography.X509Certificates;
using MathNet.Numerics.Distributions;


namespace CanvasPedigreeCaller
{
    public class CopyNumberModel
    {
        public List<int> CopyNumber = new List<int>();
        public List<Tuple<int,int>> Genotypes = new List<Tuple<int, int>>();
        public List<NegativeBinomial> CnDistribution = new List<NegativeBinomial>();
        public List<NegativeBinomial> Allele1Distribution = new List<NegativeBinomial>();
        public List<NegativeBinomial> Allele2Distribution = new List<NegativeBinomial>();

        public NegativeBinomial NegativeBinomialWrapper(double mean, double variance)
        {
            double r = Math.Pow(mean, 2)/Math.Max(variance - mean, 1.0);
            double p = r/(r + mean);
            return new NegativeBinomial(r, p);
        }

        public CopyNumberModel(int numCnStates, double haploidMean, double variance)
        {
            for (int cn = 0; cn < numCnStates; cn++)
            {
                CnDistribution.Add(NegativeBinomialWrapper(haploidMean*cn, variance));
                CopyNumber.Add(cn);
            }
        }

        public CopyNumberModel(int numCnStates, double haploidMeanAllele1, double haploidMeanAllele2, double variance)
        {
            for (int cn = 0; cn < numCnStates; cn++)
            {
                for (int gt = 0; gt <= cn; cn++)
                {
                    Allele1Distribution.Add(NegativeBinomialWrapper(haploidMeanAllele1 * gt, variance));
                    Allele2Distribution.Add(NegativeBinomialWrapper(haploidMeanAllele2 * cn-gt, variance));
                    Genotypes.Add(new Tuple<int, int>(gt, cn - gt));
                }
            }
        }

        public List<double> GetCnLikelihood(double y)
        {
            return CnDistribution.Select(x => x.Probability(Convert.ToInt32(y))).ToList();
        }

        public List<double> GetGtLikelihood(Tuple<int,int> gt)
        {
            List<double> lieklyhood = new List<double>();
            int size = Allele1Distribution.Count;
            for (int i = 0; i < size; i++)
                lieklyhood.Add(Allele1Distribution[i].Probability(gt.Item1) * Allele2Distribution[i].Probability(gt.Item2));
            return lieklyhood;
        }
    }
}
