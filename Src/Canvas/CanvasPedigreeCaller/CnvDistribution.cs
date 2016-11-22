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

        public List<double> GetMarginalProbability(int[] indices, string sampleName)
        {
            int sampleIndex = GetSampleIndex(sampleName);
            var marginalProbability = Enumerable.Repeat(0.0, indices.Length).ToList();
            double normalizationFactor = 0.0;
            var allSampleIndeces = new Permutations<int>(Enumerable.Range(0, indices.Length).ToList(), GenerateOption.WithoutRepetition);
            foreach (int copyNumberState in Enumerable.Range(0, indices.Length))
            {
                foreach (var copyNumberIndex in allSampleIndeces.Where(x => x[sampleIndex] == copyNumberState).ToArray())
                {
                    var index = (int[]) copyNumberIndex;
                    marginalProbability[copyNumberState] += GetJointProbability(index);
                    normalizationFactor += GetJointProbability(index);
                }
            }
            var query = marginalProbability.Select(x => {x = x/normalizationFactor; return x;});
            return marginalProbability;
        }

        public List<double> GetConditionalProbability(int[] indices, string sampleName, int sampleValue)
        {
            int sampleIndex = GetSampleIndex(sampleName);
            var conditionalProbability = Enumerable.Repeat(0.0, indices.Length).ToList();
            double normalizationFactor = 0.0;
            var allSampleIndeces = new Permutations<int>(Enumerable.Range(0, indices.Length).ToList(), GenerateOption.WithoutRepetition);
            foreach (int copyNumberState in Enumerable.Range(0, indices.Length))
            {
                foreach (var copyNumberIndex in allSampleIndeces.Where(x => x[sampleIndex] == copyNumberState).ToArray())
                {
                    var index = (int[])copyNumberIndex;
                    conditionalProbability[copyNumberState] += GetJointProbability(index);
                    normalizationFactor += GetJointProbability(index);
                }
            }
            var query = conditionalProbability.Select(x => { x = x / normalizationFactor; return x; });
            return conditionalProbability;
        }

    }
}
