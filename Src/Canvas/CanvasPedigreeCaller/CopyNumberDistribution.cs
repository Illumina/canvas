using System;
using System.Collections.Generic;
using System.Linq;

namespace CanvasPedigreeCaller
{
    public class CopyNumberDistribution
    {
        private readonly Array _probability;
        public List<string> SampleNames { get; }
        public int Count { get; }
        public List<int[]> Indices { get; }

        public CopyNumberDistribution(int nCopyNumbers, List<string> names)
        {
            int nSamples = names.Count;
            if (names.Count != nSamples)
                throw new ArgumentException($"List of sample names must be equal to the number of samples ");
            var dimensionSizes = Enumerable.Repeat(nCopyNumbers, nSamples).ToArray();
            _probability = Array.CreateInstance(typeof(double), dimensionSizes);
            SampleNames = names;
            Count = nSamples;
            Indices = new List<int[]>();
        }

        private int GetSampleIndex(string sampleName)
        {
            return SampleNames.FindIndex(x => x.StartsWith(sampleName));
        }

        public double GetJointProbability(int [] indices)
        {
            return Convert.ToDouble(_probability.GetValue(indices));
        }

        public void SetJointProbability(double probability, int[] indices, bool skipIndex = false)
        {
            _probability.SetValue(probability, indices);
            if (!skipIndex && !Indices.Exists(index => index.SequenceEqual(indices)) && probability > 0)
                Indices.Add(indices);
        }

        public List<double> GetMarginalProbability(int nSamples, int nCopies, string sampleName)
        {
            int sampleIndex = GetSampleIndex(sampleName);
            var marginalProbability = Enumerable.Repeat(0.0, nCopies).ToList();

            double normalizationFactor = 0.0;
            foreach (int copyNumberState in Enumerable.Range(0, nCopies))
            {
                foreach (var copyNumberIndex in Indices.Where(x => x[sampleIndex] == copyNumberState).ToArray())
                {
                    if (GetJointProbability(copyNumberIndex.ToArray()) > 0.0)
                    {
                        marginalProbability[copyNumberState] += GetJointProbability(copyNumberIndex.ToArray());
                        normalizationFactor += GetJointProbability(copyNumberIndex.ToArray());
                    }
                }
            }
            marginalProbability.Select(x => {x = x/normalizationFactor; return x;});
            return marginalProbability;
        }



        public void SetConditionalProbability(int nSamples, int nCopies, string sampleName, int sampleCnValue, double sampleMarginalProbability)
        {
            int sampleIndex = GetSampleIndex(sampleName);

            foreach (var index in Indices.Where(x => x[sampleIndex] == sampleCnValue))
            {
                SetJointProbability(GetJointProbability(index.ToArray()) / sampleMarginalProbability, index.ToArray(), true);
            }
        }
    }
}
