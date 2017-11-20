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
        public double MaximalLikelihood;
        private readonly Dictionary<ISampleMap<Genotype>, double> _jointLikelihoods;
        private double _totalLikelihood;

        public JointLikelihoods()
        {
            MaximalLikelihood = double.MinValue;
            _jointLikelihoods = new Dictionary<ISampleMap<Genotype>, double>();
            _totalLikelihood = 0;
        }

        public void SetJointLikelihood(ISampleMap<Genotype> samplesGenotypes, double likelihood)
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
    }
}