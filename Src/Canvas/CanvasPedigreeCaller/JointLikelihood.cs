using System.Collections.Generic;
using System.IO;
using System.Linq;
using CanvasCommon;
using Illumina.Common.FileSystem;
using Isas.Framework.DataTypes;
using Isas.Framework.DataTypes.Maps;

namespace CanvasPedigreeCaller
{
    internal class JointLogLikelihoods
    {
        public double MaximalLogLikelihood;
        private readonly Dictionary<ISampleMap<Genotype>, double> _jointLogLikelihoods;
        private double _totalLogLikelihood;

        public JointLogLikelihoods()
        {
            MaximalLogLikelihood = double.MinValue;
            var comparer = new SampleGenotypeComparer();
            _jointLogLikelihoods = new Dictionary<ISampleMap<Genotype>, double>(comparer);
            _totalLogLikelihood = 0;
        }

        public void AddJointLikelihood(ISampleMap<Genotype> samplesGenotypes, double likelihood)
        {
            _jointLogLikelihoods[samplesGenotypes] = likelihood;
            _totalLogLikelihood += likelihood;
        }

        public double GetJointLikelihood(ISampleMap<Genotype> samplesGenotypes)
        {
            return _jointLogLikelihoods[samplesGenotypes];
        }

        public double GetMarginalLikelihood(KeyValuePair<SampleId, Genotype> samplesGenotype)
        {
            return _jointLogLikelihoods.Where(kvp => Equals(kvp.Key[samplesGenotype.Key], samplesGenotype.Value)).Select(kvp => kvp.Value).Sum() / _totalLogLikelihood;
        }

        public double GetMarginalNonAltLikelihood(KeyValuePair<SampleId, Genotype> samplesGenotype)
        {
            return _jointLogLikelihoods.Where(kvp => !Equals(kvp.Key[samplesGenotype.Key], samplesGenotype.Value)).Select(kvp => kvp.Value).Sum() / _totalLogLikelihood;
        }

        private class SampleGenotypeComparer : IEqualityComparer<ISampleMap<Genotype>>
        {
            public bool Equals(ISampleMap<Genotype> x, ISampleMap<Genotype> y)
            {
                return x.SequenceEqual(y);
            }

            public int GetHashCode(ISampleMap<Genotype> obj)
            {
                return obj.Aggregate(17, (hash, value) => hash + value.GetHashCode() * 31);
            }
        }
    }
}