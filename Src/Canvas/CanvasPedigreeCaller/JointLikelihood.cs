using System.Collections.Generic;
using System.IO;
using System.Linq;
using CanvasCommon;
using Illumina.Common.FileSystem;
using Isas.Framework.DataTypes;
using Isas.Framework.DataTypes.Maps;

namespace CanvasPedigreeCaller
{
    internal class JointLikelihoods
    {
        public double MaximalLogLikelihood;
        private readonly Dictionary<ISampleMap<Genotype>, double> _jointLikelihoods;
        private double _totalMarginalLikelihood;

        public JointLikelihoods()
        {
            MaximalLogLikelihood = double.NegativeInfinity;
            var comparer = new SampleGenotypeComparer();
            _jointLikelihoods = new Dictionary<ISampleMap<Genotype>, double>(comparer);
            _totalMarginalLikelihood = 0;
        }

        public void AddJointLikelihood(ISampleMap<Genotype> samplesGenotypes, double likelihood)
        {
            _jointLikelihoods[samplesGenotypes] = likelihood;
            _totalMarginalLikelihood += likelihood;
        }

        public double GetJointLikelihood(ISampleMap<Genotype> samplesGenotypes)
        {
            return _jointLikelihoods[samplesGenotypes];
        }

        // in a pedigree with the map (SampleId[M]=>Genotype[G], M: parents, offspring, G: genotype), estimate posterior likelihood as
        // (SampleId[M]=>Genotype[G], sum over all M=m, G=g)/(SampleId[M]=>Genotype[G], sum over all M and G, i.e. what is the probability of 
        // pedigree member X having genotype Y
        public double GetMarginalLikelihood(KeyValuePair<SampleId, Genotype> samplesGenotype)
        {
            return _jointLikelihoods.Where(kvp => Equals(kvp.Key[samplesGenotype.Key], samplesGenotype.Value)).Select(kvp => kvp.Value).Sum() /
                _totalMarginalLikelihood;
        }

        // in a pedigree with the map (SampleId[M]=>Genotype[G], M: parents, offspring, G: genotype), estimate posterior likelihood as
        // (SampleId[M]=>Genotype[G], sum over all M!=m, G!=g)/(SampleId[M]=>Genotype[G], sum over all M and G, i.e. what is the probability of 
        // pedigree member X not having genotype Y
        public double GetMarginalNonAltLikelihood(KeyValuePair<SampleId, Genotype> samplesGenotype)
        {
            return _jointLikelihoods.Where(kvp => !Equals(kvp.Key[samplesGenotype.Key], samplesGenotype.Value)).Select(kvp => kvp.Value).Sum() / 
                _totalMarginalLikelihood;
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