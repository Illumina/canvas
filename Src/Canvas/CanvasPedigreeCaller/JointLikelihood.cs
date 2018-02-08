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
        public double MaximalLikelihood;
        private readonly Dictionary<ISampleMap<Genotype>, double> _jointLikelihoods;
        private double _totalLikelihood;

        public JointLogLikelihoods()
        {
            MaximalLikelihood = double.MinValue;
            var comparer = new SampleGenotypeComparer();
            _jointLikelihoods = new Dictionary<ISampleMap<Genotype>, double>(comparer);
            _totalLikelihood = 0;
        }

        public void AddJointLikelihood(ISampleMap<Genotype> samplesGenotypes, double likelihood)
        {
            _jointLikelihoods[samplesGenotypes] = likelihood;
            _totalLikelihood += likelihood;
        }

        public double GetJointLikelihood(ISampleMap<Genotype> samplesGenotypes)
        {
            return _jointLikelihoods[samplesGenotypes];
        }

        public double GetMarginalLikelihood(KeyValuePair<SampleId, Genotype> samplesGenotype)
        {
            return _jointLikelihoods.Where(kvp => Equals(kvp.Key[samplesGenotype.Key], samplesGenotype.Value)).Select(kvp => kvp.Value).Sum() / _totalLikelihood;
        }

        public double GetMarginalNonAltLikelihood(KeyValuePair<SampleId, Genotype> samplesGenotype)
        {
            return _jointLikelihoods.Where(kvp => !Equals(kvp.Key[samplesGenotype.Key], samplesGenotype.Value)).Select(kvp => kvp.Value).Sum() / _totalLikelihood;
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