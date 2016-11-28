using System;
using System.Collections.Generic;
using System.Dynamic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Combinatorics.Collections;
using MathNet.Numerics;

namespace CanvasPedigreeCaller
{
    public class CopyNumberDistribution
    {
        private readonly Array _probability;
        public List<string> SampleNames { get; }
        public int Count { get; }


        public CopyNumberDistribution(int nCopyNumbers, List<string> names)
        {
            int nSamples = names.Count;
            if (names.Count != nSamples)
                throw new ArgumentException($"List of sample names must be equal to the number of samples ");
            var dimensionSizes = Enumerable.Repeat(nCopyNumbers, nSamples).ToArray();
            _probability = Array.CreateInstance(typeof(double), dimensionSizes);
            SampleNames = names;
            Count = nSamples;
        }

        private int GetSampleIndex(string sampleName)
        {
            return SampleNames.FindIndex(x => x.StartsWith(sampleName));
        }

        public double GetJointProbability(int [] indices)
        {
            return Convert.ToDouble(_probability.GetValue(indices));
        }

        public void SetJointProbability(double probability, int[] indices)
        {
            _probability.SetValue(probability, indices);
        }

        public List<double> GetMarginalProbability(int nSamples, int nCopies, string sampleName)
        {
            int sampleIndex = GetSampleIndex(sampleName);
            var marginalProbability = Enumerable.Repeat(0.0, nCopies).ToList();

            double normalizationFactor = 0.0;
            var copies = Enumerable.Range(0, nCopies).ToList();
            var combinations = new Combinations<int>(copies, nSamples, GenerateOption.WithRepetition);
            var allSampleIndeces = new List<List<int>>();

            foreach (var combination in combinations)
            {
                var permutation = new Permutations<int>(combination, GenerateOption.WithoutRepetition);
                allSampleIndeces.AddRange(permutation.Select(permutationElement => permutationElement.ToList()));
            }

            foreach (int copyNumberState in Enumerable.Range(0, nCopies))
            {
                foreach (var copyNumberIndex in allSampleIndeces.Where(x => x[sampleIndex] == copyNumberState).ToArray())
                {
                    if (GetJointProbability(copyNumberIndex.ToArray()) > 0.0)
                    {
                        marginalProbability[copyNumberState] += GetJointProbability(copyNumberIndex.ToArray());
                        normalizationFactor += GetJointProbability(copyNumberIndex.ToArray());
                    }
                }
            }
            var query = marginalProbability.Select(x => {x = x/normalizationFactor; return x;});
            return marginalProbability;
        }



        public void SetConditionalProbability(int nSamples, int nCopies, string sampleName, int sampleCnValue, double sampleMarginalProbability)
        {
            int sampleIndex = GetSampleIndex(sampleName);
            var copies = Enumerable.Range(0, nCopies).ToList();
            var combinations = new Combinations<int>(copies, nSamples, GenerateOption.WithRepetition);
            var allSampleIndeces = new List<List<int>>();

            foreach (var combination in combinations)
            {
                var permutation = new Permutations<int>(combination, GenerateOption.WithoutRepetition);
                allSampleIndeces.AddRange(permutation.Select(permutationElement => permutationElement.ToList()));
            }
            foreach (var index in allSampleIndeces.Where(x => x[sampleIndex] == sampleCnValue))
            {
                SetJointProbability(GetJointProbability(index.ToArray()) / sampleMarginalProbability, index.ToArray());
            }
        }
    }
}
