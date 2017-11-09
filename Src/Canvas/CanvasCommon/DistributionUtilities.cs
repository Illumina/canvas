using System;
using System.Collections.Generic;
using System.Linq;
using Combinatorics.Collections;
using MathNet.Numerics;

namespace CanvasCommon
{
    public static class DistributionUtilities
    {
        public static List<List<int>> GetGenotypeCombinations(int numberOfStates, int currentState)
        {
            const int diploidState = 2;
            const int maxNumberOfStates = 4;
            if (numberOfStates > maxNumberOfStates)
            {
                Console.WriteLine($"CanvasPartition.exe: number of states/samples {numberOfStates} is too large for full combinatorial enumeration\n" +
                    $"Settings to the default number {maxNumberOfStates}");
                numberOfStates = maxNumberOfStates;
            }

            double upperSetBound = SpecialFunctions.Factorial(numberOfStates) * SpecialFunctions.Factorial(numberOfStates / 2);
            var allCombinations = new List<List<int>>(Convert.ToInt32(upperSetBound));
            for (int numberOfDiploidStates = 1; numberOfDiploidStates < numberOfStates; numberOfDiploidStates++)
            {
                var states = Enumerable.Repeat(currentState, numberOfStates - numberOfDiploidStates)
                    .Concat(Enumerable.Repeat(diploidState, numberOfDiploidStates));
                var permutations = new Permutations<int>(states.ToList(), GenerateOption.WithoutRepetition);
                var list = permutations.Select(x => x.ToList()).ToList();
                allCombinations.AddRange(list);
            }
            return allCombinations.Count == 0 ? new List<List<int>>{new List<int>{currentState}} : allCombinations;
        }

        /// <summary>
        /// Calculate density using the negative binomial distribution. clumpingParameter reflects the tightness of 
        /// distibution.
        /// </summary>
        /// <param name="mean"></param>
        /// <param name="variance"></param>
        /// <param name="maxValue"></param>
        /// <param name="adjustClumpingParameter"></param>
        /// <returns></returns>
        public static List<double> NegativeBinomialWrapper(double mean, double variance, int maxValue, bool adjustClumpingParameter = false)
        {
            var negativeBinomialDensity = Enumerable.Repeat(0.0, maxValue).ToList();
            const double minMeanValue = 0.1; // handle the case of zero counts
            const double adjustClumpingParameterMin1 = 6.0; 
            const double adjustClumpingParameterMin2 = 2.0; 

            // handle the case where variance is less than mean
            double clumpingParameter = Math.Pow(Math.Max(mean, minMeanValue), 2) / (Math.Max(variance, mean * 1.2) - mean);
            clumpingParameter = Math.Max(adjustClumpingParameter ? adjustClumpingParameterMin1 : adjustClumpingParameterMin2, clumpingParameter);
            for (int x = 0; x < maxValue; x++)
            {
                // density function for negative binomial distribution
                var tmpDensity = Math.Exp(Math.Log(Math.Pow(1 + mean / clumpingParameter, -clumpingParameter)) + Math.Log(Math.Pow(mean / (mean + clumpingParameter), x)) + SpecialFunctions.GammaLn(clumpingParameter + x) -
                             SpecialFunctions.FactorialLn(x) - SpecialFunctions.GammaLn(clumpingParameter));
                negativeBinomialDensity[x] = Double.IsNaN(tmpDensity) || Double.IsInfinity(tmpDensity) ? 0 : tmpDensity;
            }
            return negativeBinomialDensity;
        }
    }
}
