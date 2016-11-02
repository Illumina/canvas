using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography.X509Certificates;
using CanvasCommon;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;


namespace CanvasPedigreeCaller
{
    public class CopyNumberModel
    {
        public List<Tuple<int,int>> Genotypes = new List<Tuple<int, int>>();
        public List<NegativeBinomial> CnDistribution = new List<NegativeBinomial>();
        Tuple<NegativeBinomial, NegativeBinomial> [][] AlleleDistribution;

        public NegativeBinomial NegativeBinomialWrapper(double mean, double variance)
        {
            double r = Math.Pow(Math.Max(mean, 0.1), 2)/Math.Max(variance - mean, 1.0);
            double p = r/(r + mean);
            return new NegativeBinomial(r, p);
        }

        public CopyNumberModel(int numCnStates, double haploidMean, double haploidMafMean, double variance, double mafVariance)
        {
            for (int cn = 0; cn < numCnStates; cn++)
            {
                CnDistribution.Add(NegativeBinomialWrapper(haploidMean*cn, variance));
            }
            AlleleDistribution = new Tuple<NegativeBinomial, NegativeBinomial>[numCnStates][];
            for (int i = 0; i < numCnStates; i++)
                AlleleDistribution[i] = new Tuple<NegativeBinomial, NegativeBinomial>[numCnStates];

            for (int gt1 = 0; gt1 < numCnStates; gt1++)
            {
                for (int gt2 = 0; gt2 < numCnStates; gt2++)
                {
                    var gt1Prob = NegativeBinomialWrapper(haploidMafMean * gt1, mafVariance);
                    var gt2Prob = NegativeBinomialWrapper(haploidMafMean * gt2, mafVariance);
                    AlleleDistribution[gt1][gt2] = new Tuple<NegativeBinomial, NegativeBinomial>(gt1Prob, gt2Prob);
                }
            }
        }


        public List<double> GetCnLikelihood(double y)
        {
            return CnDistribution.Select(x => x.Probability(Convert.ToInt32(y))).ToList();
        }

        public double[][] GetGtLikelihood(Tuple<int,int> gt)
        {
            int nrows = AlleleDistribution.Length;
            int ncols = AlleleDistribution.First().Length;
            double[][] lieklyhood = Utilities.MatrixCreate(nrows, ncols);

            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    if (AlleleDistribution[i][j] != null)
                        lieklyhood[i][j] = AlleleDistribution[i][j].Item1.Probability(gt.Item1)*
                                           AlleleDistribution[i][j].Item2.Probability(gt.Item2);
                    else
                        lieklyhood[i][j] = 0;
                }
            }
            return lieklyhood;
        }
    }
}
