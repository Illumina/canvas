using System;
using System.Collections.Generic;
using System.Linq;
using Isas.Framework.DataTypes;

namespace CanvasPedigreeCaller
{
    public class CopyNumbersLikelihoods
    {
        private readonly Array _probability;
        public List<string> SampleNames { get; }
        public int Count { get; }
        public List<int[]> Indices { get; }
        public double MaximalLikelihood { get; set; }
        public SampleList<Dictionary<int, double>> SingleSampleLikelihoods { get; }

        public CopyNumbersLikelihoods(SampleList<Dictionary<int, double>> singleSampleLikelihoods, int nCopyNumbers)
        {
            this.SingleSampleLikelihoods = singleSampleLikelihoods;
            Indices = new List<int[]>();
            this.SampleNames = singleSampleLikelihoods.SampleIds.Select(id => id.ToString()).ToList();
            var dimensionSizes = Enumerable.Repeat(nCopyNumbers, this.SampleNames.Count).ToArray();
            _probability = Array.CreateInstance(typeof(double), dimensionSizes);
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

        public List<double> GetMarginalProbability(int nCopies, string sampleName)
        {
            int sampleIndex = GetSampleIndex(sampleName);
            var marginalProbability = Enumerable.Repeat(0.0, nCopies).ToList();

            foreach (int copyNumberState in Enumerable.Range(0, nCopies))
            {
                foreach (var copyNumberIndex in Indices.Where(x => x[sampleIndex] == copyNumberState).ToArray())
                {
                    if (GetJointProbability(copyNumberIndex.ToArray()) > 0.0)
                        marginalProbability[copyNumberState] += GetJointProbability(copyNumberIndex.ToArray());
                }
            }
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
