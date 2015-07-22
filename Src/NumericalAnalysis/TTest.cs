using System;

namespace Illumina.NumericalAnalysis
{
    public class Test
    {
        /// <summary>
        ///     Returns Sp^2
        /// </summary>
        /// <param name="n1">num samples</param>
        /// <param name="n2"></param>
        /// <param name="s12">sigma^2</param>
        /// <param name="s22"></param>
        /// <returns></returns>
        public static double PooledEstimatorForSigma(double n1, double n2, double s12, double s22)
        {
            return (((n1 - 1)*s12) + ((n2 - 1)*s22))/(n1 + n2 - 2);
        }

        public static double TwoPopulationTTest(double m1, double m2,
                                                double n1, double n2, double sp)
        {
            return (m1 - m2)/(sp*Math.Sqrt((1/n1) + (1/n2)));
        }

        /// <summary>
        /// </summary>
        /// <param name="m1">Support</param>
        /// <param name="m2"></param>
        /// <param name="n1">Coverage</param>
        /// <param name="n2"></param>
        /// <param name="levelOfSignificance"> </param>
        /// <returns></returns>
        public static double[] GetTValue(double m1, double m2, double n1, double n2, double levelOfSignificance)
        {
            double degreesOfFreedom = n1 + n2 - 2;

            //if we have  very low depth, refuse to do anything.
            //If we do need to go lower, then i would need to input values for the lookup table here.
            if (degreesOfFreedom < 30)
            {
                return new[] {double.NaN, double.NaN};
            }

            double sp2 = PooledEstimatorForSigma(n1, n2, m1, m2);


            double tvalue = TwoPopulationTTest(m1, m2, n1, n2, sp2);

            return new[] {tvalue, degreesOfFreedom};
        }

        private static bool ValueAcceptable(double levelOfSignificance, double tvalue, double degreesOfFreedom)
        {
            double alphaOver2 = levelOfSignificance/2.0;
            double rejectionRegion = 1.282;

            if (degreesOfFreedom < 30)
            {
                //throw new Exception("Its best if you go look this up in a table...");
                return false; //just don't call anything for now.
            }

            //From "Mathematical Statistics With Applications" 6th ed.
            if (alphaOver2 < 0.005)
                rejectionRegion = 2.576;
            else if (alphaOver2 < 0.01)
                rejectionRegion = 2.576;
            else if (alphaOver2 < 0.025)
                rejectionRegion = 2.326;
            else if (alphaOver2 < 0.05)
                rejectionRegion = 1.960;
            else if (alphaOver2 < 0.1)
                rejectionRegion = 1.645;

            if (Math.Abs(tvalue) > rejectionRegion)
            {
                return false;
            }

            return true;
        }
    }
}