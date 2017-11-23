using System;
using System.Collections.Generic;
using System.Linq;
using CanvasCommon;
using Isas.Framework.DataTypes;
using Isas.Framework.DataTypes.Maps;
using MathNet.Numerics.Distributions;

namespace CanvasPedigreeCaller
{
    class PedigreeInfo
    {
        public List<SampleId> ParentsIds { get; }
        public List<SampleId> OffspringIds { get; }
        public List<List<Genotype>> OffspringPhasedGenotypes { get; }
        public List<List<Genotype>> OffspringTotalCopyNumberGenotypes { get; }

        public double[][] TransitionMatrix { get; }

        private PedigreeInfo(List<SampleId> offspringIds, List<SampleId> parentsIds, List<List<Genotype>> offspringPhasedGenotypes,
            List<List<Genotype>> offspringTotalCopyNumberGenotypes, double[][] transitionMatrix)
        {
            OffspringIds = offspringIds;
            ParentsIds = parentsIds;
            OffspringPhasedGenotypes = offspringPhasedGenotypes;
            OffspringTotalCopyNumberGenotypes = offspringTotalCopyNumberGenotypes;
            TransitionMatrix = transitionMatrix;
        }

        public static PedigreeInfo GetPedigreeInfo(ISampleMap<SampleType> kinships, PedigreeCallerParameters callerParameters)
        {
            var parentsIds = kinships.Where(kin => kin.Value.Equals(SampleType.Father) || kin.Value.Equals(SampleType.Mother)).Select(kin => kin.Key)
                .ToList();
            var offspringIds = kinships.Where(kin => kin.Value.Equals(SampleType.Proband) || kin.Value.Equals(SampleType.Sibling)).Select(kin => kin.Key)
                .ToList();
            var parentalPhasedGenotypes = GeneratePhasedGenotype(callerParameters.MaximumCopyNumber);
            var parentalTotalCopyNumberGenotypes = Enumerable.Range(0, callerParameters.MaximumCopyNumber).Select(Genotype.Create).ToList();
            var offspringPhasedGenotypes = GetOffspringGenotypes(callerParameters, parentalPhasedGenotypes, offspringIds);
            var offspringTotalCopyNumberGenotypes = GetOffspringGenotypes(callerParameters, parentalTotalCopyNumberGenotypes, offspringIds);
            var transitionMatrix = GetTransitionMatrix(callerParameters.MaximumCopyNumber);
            return new PedigreeInfo(offspringIds, parentsIds, offspringPhasedGenotypes, offspringTotalCopyNumberGenotypes, transitionMatrix);
        }

        private static List<List<Genotype>> GetOffspringGenotypes(PedigreeCallerParameters callerParameters, List<Genotype> genotypes, List<SampleId> offspringIds)
        {
            var offspringGenotypes = new List<List<Genotype>>(Convert.ToInt32(Math.Pow(genotypes.Count, offspringIds.Count)));
            GenerateOffspringGenotypes(offspringGenotypes, genotypes, offspringIds.Count, new List<Genotype>());
            if (offspringGenotypes.Count > callerParameters.MaxNumOffspringGenotypes)
            {
                offspringGenotypes.Shuffle();
                offspringGenotypes = offspringGenotypes.Take(callerParameters.MaxNumOffspringGenotypes).ToList();
            }
            return offspringGenotypes;
        }

        public static List<Genotype> GeneratePhasedGenotype(int numCnStates)
        {
            var genotypes = new List<Genotype>();
            for (int cn = 0; cn < numCnStates; cn++)
            {
                for (int gt = 0; gt <= cn; gt++)
                {
                    genotypes.Add(Genotype.Create(new PhasedGenotype(gt, cn - gt)));
                }
            }
            return genotypes;
        }

        public static void GenerateOffspringGenotypes(List<List<Genotype>> offspringGenotypes, List<Genotype> genotypeSet, int nOffsprings, List<Genotype> partialGenotypes)
        {
            if (nOffsprings > 0)
            {
                foreach (var genotype in genotypeSet)
                {
                    GenerateOffspringGenotypes(offspringGenotypes, genotypeSet, nOffsprings - 1, 
                        partialGenotypes.Concat(new List<Genotype> { genotype }).ToList());
                }
            }
            if (nOffsprings == 0)
            {
                offspringGenotypes.Add(partialGenotypes);
            }
        }

        public static double[][] GetTransitionMatrix(int numCnStates)
        {
            double[][] transitionMatrix = Utilities.MatrixCreate(numCnStates, numCnStates);
            transitionMatrix[0][0] = 1.0;
            for (int gt = 1; gt < numCnStates; gt++)
                transitionMatrix[0][gt] = 0;

            for (int cn = 1; cn < numCnStates; cn++)
            {
                var gtLikelihood = new Poisson(Math.Max(cn / 2.0, 0.1));
                for (int gt = 0; gt < numCnStates; gt++)
                    transitionMatrix[cn][gt] = gtLikelihood.Probability(gt);
            }
            return transitionMatrix;
        }
    }
}
