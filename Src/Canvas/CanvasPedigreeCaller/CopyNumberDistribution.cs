using System;
using System.Collections.Generic;
using System.Linq;
using Isas.Framework.DataTypes;

namespace CanvasPedigreeCaller
{
    public class CopyNumbersLikelihoods
    {
        private readonly Array _jointLikelihoods;
        public List<SampleList<int>> CopyNumbers { get; }
        public double MaximalLikelihood { get; set; }
        public SampleList<Dictionary<int, double>> SingleSampleLikelihoods { get; }

        public CopyNumbersLikelihoods(SampleList<Dictionary<int, double>> singleSampleLikelihoods, int nCopyNumbers)
        {
            SingleSampleLikelihoods = singleSampleLikelihoods;
            CopyNumbers = new List<SampleList<int>>();
            var dimensionSizes = Enumerable.Repeat(nCopyNumbers, singleSampleLikelihoods.Count()).ToArray();
            _jointLikelihoods = Array.CreateInstance(typeof(double), dimensionSizes);
        }

        public double GetJointLikelihood(SampleList<int> copyNumbers)
        {
            var copyNumberIndices = copyNumbers.SampleData.ToArray();
            return Convert.ToDouble(_jointLikelihoods.GetValue(copyNumberIndices));
        }

        public void SetJointLikelihood(double probability, SampleList<int> copyNumbers, bool skipIndex = false)
        {
            var copyNumberIndices = copyNumbers.SampleData.ToArray();
            _jointLikelihoods.SetValue(probability, copyNumberIndices);
            if (!skipIndex && !CopyNumbers.Exists(index => index.SequenceEqual(copyNumbers)) && probability > 0)
                CopyNumbers.Add(copyNumbers);
        }

        public List<double> GetMarginalProbability(int nCopies, SampleId sampleId)
        {
            var marginalProbability = Enumerable.Repeat(0.0, nCopies).ToList();

            foreach (int copyNumberState in Enumerable.Range(0, nCopies))
            {
                foreach (var copyNumberIndex in CopyNumbers.Where(x => x[sampleId] == copyNumberState).ToArray())

                    marginalProbability[copyNumberState] += GetJointLikelihood(copyNumberIndex);
            }
            return marginalProbability;
        }

    }
}
