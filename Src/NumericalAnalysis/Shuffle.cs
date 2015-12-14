using System;
using System.Collections.Generic;

namespace Illumina.NumericalAnalysis
{
    public static class Shuffle
    {
        public static void FisherYates<T>(List<T> list)
        {
            Random rnd = new Random();

            for (int i = list.Count; i > 1; i--)
            {
                // Pick random element to swap.
                int j = rnd.Next(i);

                // swap the elements
                T tmp = list[j];
                list[j] = list[i - 1];
                list[i - 1] = tmp;
            }
        }
    }
}