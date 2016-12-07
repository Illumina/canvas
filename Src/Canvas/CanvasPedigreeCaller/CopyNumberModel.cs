using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography.X509Certificates;
using CanvasCommon;
using MathNet.Numerics;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;


namespace CanvasPedigreeCaller
{
    public class CopyNumberModel
    {
        public List<Tuple<int,int>> Genotypes = new List<Tuple<int, int>>();
        public List<List<double>> CnDistribution = new List<List<double>>();
        readonly Tuple<List<double>, List<double>> [][] _alleleDistribution;

        public static List<double> NegativeBinomialWrapper(double mean, double variance, int maxValue)
        {
            var density = Enumerable.Repeat(0.0, maxValue).ToList();
            double r = Math.Pow(Math.Max(mean, 0.1), 2) / (Math.Max(variance, mean * 1.2) - mean);
            for (int x = 0; x < maxValue; x++)
            {
                var tmpDensity = Math.Exp(Math.Log(Math.Pow(1 + mean / r, -r)) + Math.Log(Math.Pow(mean / (mean + r), x)) + SpecialFunctions.GammaLn(r + x) -
                             SpecialFunctions.FactorialLn(x) - SpecialFunctions.GammaLn(r));
                density[x] = Double.IsNaN(tmpDensity) || Double.IsInfinity(tmpDensity) ? 0 : tmpDensity;
            }
            return density;
        }
        public CopyNumberModel(int numCnStates, double haploidMean, double haploidMafMean, double variance, double mafVariance, int maxValue)
        {
            for (int cn = 0; cn < numCnStates; cn++)
            {
                var multiplier = 1.0;
                if (cn == 1)
                    multiplier = 0.9;
                if (cn == 3)
                    multiplier = 1.1;

                CnDistribution.Add(NegativeBinomialWrapper(haploidMean*cn*multiplier, variance, maxValue));
            }

            _alleleDistribution = new Tuple<List<double>, List<double>>[numCnStates][];

            for (int i = 0; i < numCnStates; i++)
                _alleleDistribution[i] = new Tuple<List<double>, List<double>>[numCnStates];

            for (int gt1 = 0; gt1 < numCnStates; gt1++)
            {
                for (int gt2 = 0; gt2 < numCnStates; gt2++)
                {
                    var gt1Probabilities = NegativeBinomialWrapper(haploidMafMean * gt1, mafVariance, maxValue);
                    var gt2Probabilities = NegativeBinomialWrapper(haploidMafMean * gt2, mafVariance, maxValue);
                    _alleleDistribution[gt1][gt2] = new Tuple<List<double>, List<double>>(gt1Probabilities, gt2Probabilities);
                }
            }
        }


        public List<double> GetCnLikelihood(double y)
        {
            return CnDistribution.Select(x => x[Convert.ToInt32(y)]).ToList();
        }

        public double[][] GetGtLikelihood(Tuple<int,int> gt)
        {
            int nrows = _alleleDistribution.Length;
            int ncols = _alleleDistribution.First().Length;
            double[][] lieklyhood = Utilities.MatrixCreate(nrows, ncols);

            for (int i = 0; i < nrows; i++)
            {
                for (int j = 0; j < ncols; j++)
                {
                    if (_alleleDistribution[i][j] != null)
                        lieklyhood[i][j] = _alleleDistribution[i][j].Item1[gt.Item1]*
                                           _alleleDistribution[i][j].Item2[gt.Item2];
                    else
                        lieklyhood[i][j] = 0;
                }
            }
            return lieklyhood;
        }
    }
}
