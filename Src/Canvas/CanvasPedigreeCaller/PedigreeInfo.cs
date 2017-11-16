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
        public List<List<PhasedGenotype>> OffspringGenotypes { get; }
        public double[][] TransitionMatrix { get; }

        private PedigreeInfo(List<SampleId> offspringIds, List<SampleId> parentsIds, List<List<PhasedGenotype>> offspringGenotypes, double[][] transitionMatrix)
        {
            OffspringIds = offspringIds;
            ParentsIds = parentsIds;
            OffspringGenotypes = offspringGenotypes;
            TransitionMatrix = transitionMatrix;
        }

        public static PedigreeInfo GetPedigreeInfo(ISampleMap<SampleType> kinships, PedigreeCallerParameters callerParameters)
        {
            var parentsIds = kinships.Where(kin => kin.Value.Equals(SampleType.Father) || kin.Value.Equals(SampleType.Mother)).Select(kin => kin.Key)
                .ToList();
            var offspringIds = kinships.Where(kin => kin.Value.Equals(SampleType.Proband) || kin.Value.Equals(SampleType.Sibling)).Select(kin => kin.Key)
                .ToList();
            var parentalGenotypes = GenerateParentalGenotypes(callerParameters.MaximumCopyNumber);
            var offspringsGenotypes =
                new List<List<PhasedGenotype>>(Convert.ToInt32(Math.Pow(parentalGenotypes.Count, offspringIds.Count)));
            GenerateOffspringGenotypes(offspringsGenotypes, parentalGenotypes, offspringIds.Count, new List<PhasedGenotype>());
            if (offspringsGenotypes.Count > callerParameters.MaxNumOffspringGenotypes)
            {
                offspringsGenotypes.Shuffle();
                offspringsGenotypes = offspringsGenotypes.Take(callerParameters.MaxNumOffspringGenotypes).ToList();
            }
            var transitionMatrix = GetTransitionMatrix(callerParameters.MaximumCopyNumber);
            return new PedigreeInfo(offspringIds, parentsIds, offspringsGenotypes, transitionMatrix);
        }


        public static List<PhasedGenotype> GenerateParentalGenotypes(int numCnStates)
        {
            var genotypes = new List<PhasedGenotype>();
            for (int cn = 0; cn < numCnStates; cn++)
            {
                for (int gt = 0; gt <= cn; gt++)
                {
                    genotypes.Add(new PhasedGenotype(gt, cn - gt));
                }
            }
            return genotypes;
        }

        public static void GenerateOffspringGenotypes(List<List<PhasedGenotype>> offspringGenotypes, List<PhasedGenotype> genotypeSet,
            int nOffsprings, List<PhasedGenotype> partialGenotypes)
        {
            if (nOffsprings > 0)
            {
                foreach (var genotype in genotypeSet)
                {
                    GenerateOffspringGenotypes(offspringGenotypes, genotypeSet, nOffsprings - 1,
                        partialGenotypes.Concat(new List<PhasedGenotype> { genotype }).ToList());
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
