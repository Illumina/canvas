using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Isas.SequencingFiles;
using MathNet.Numerics.Distributions;


namespace CanvasPedigreeCaller
{
    public class CopyNumberModel
    {
        public List<int> CopyNumber = new List<int>();
        public List<NegativeBinomial> Distribution = new List<NegativeBinomial>();

        public NegativeBinomial NegativeBinomialWrapper(double mean, double variance)
        {
            double r = Math.Pow(mean, 2)/(Math.Pow(variance, 2) - mean);
            double p = r/(r + mean);
            return new NegativeBinomial(r, p);
        }

        public CopyNumberModel(int numCnStates, double haploidMean, double variance)
        {
            for (int cn = 0; cn < numCnStates; cn++)
            {
                Distribution.Add(NegativeBinomialWrapper(haploidMean*cn, variance));
                CopyNumber.Add(cn);
            }
        }

        public List<double> GetCnLikelihood(double y)
        {
            return Distribution.Select(x => x.Probability(Convert.ToInt32(y))).ToList();
        }

    }
}
