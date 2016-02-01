using System;

namespace Illumina.NumericalAnalysis
{
    public class Quality
    {
        public static double QtoP(double q)
        {
            return Math.Pow(10.0, -1 * q / 10f);
        }

        public static double PtoQ(double p)
        {
            return (-10.0 * Math.Log10(p));
        }
    }
}