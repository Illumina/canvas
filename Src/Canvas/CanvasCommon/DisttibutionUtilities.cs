using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.Distributions;
using System.Text;
using System.Threading.Tasks;
using Combinatorics.Collections;
using MathNet.Numerics;

namespace CanvasCommon
{
    public static class DistibutionUtilities
    {
        public static List<List<int>> GetGenotypeCombinations(int numberOfStates, int currentState)
        {
            const int diploidState = 2;
            const int maxnumberOfStates = 4;
            if (numberOfStates > maxnumberOfStates)
            {
                Console.WriteLine($"CanvasPartition.exe: number of states/samples {numberOfStates} is too large for full combinatorial enumeration\n" +
                    $"Stettign to the default number {maxnumberOfStates}");
                numberOfStates = maxnumberOfStates;
            }

            var upperSetBound = SpecialFunctions.Factorial(numberOfStates) * SpecialFunctions.Factorial(numberOfStates / 2);
            var allCombinations = new List<List<int>>(Convert.ToInt32(upperSetBound));
            for (int numberOfDiploidStates = 1; numberOfDiploidStates < numberOfStates; numberOfDiploidStates++)
            {
                var states = Enumerable.Repeat(currentState, numberOfStates - numberOfDiploidStates)
                    .Concat(Enumerable.Repeat(diploidState, numberOfDiploidStates));
                var permutations = new Permutations<int>(states.ToList(), GenerateOption.WithoutRepetition);
                var list = permutations.Select(x => x.ToList()).ToList();
                allCombinations.AddRange(list);
            }
            return allCombinations;
        }

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
    }
}
