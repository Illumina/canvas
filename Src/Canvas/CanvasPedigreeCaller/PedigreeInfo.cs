using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CanvasCommon;
using Isas.Framework.DataTypes;
using MathNet.Numerics.Distributions;

namespace CanvasPedigreeCaller
{
    class PedigreeInfo
    {
        public List<SampleId> ParentsIds { get; }
        public List<SampleId> OffspringsIds { get; }

        public List<Genotype> ParentalGenotypes { get; }
        public List<List<Genotype>> OffspringsGenotypes { get; }
        public double[][] TransitionMatrix { get; }

        public PedigreeInfo(SampleList<CanvasPedigreeCaller.Kinship> kinships, PedigreeCallerParameters callerParameters)
        {
            ParentsIds = kinships.Where(kin => kin.Value.Equals(CanvasPedigreeCaller.Kinship.Parent)).Select(kin => kin.Key)
                .ToList();
            OffspringsIds = kinships.Where(kin => kin.Value.Equals(CanvasPedigreeCaller.Kinship.Proband)).Select(kin => kin.Key)
                .ToList();
            ParentalGenotypes = GenerateParentalGenotypes(callerParameters.MaximumCopyNumber);
            OffspringsGenotypes =
                new List<List<Genotype>>(Convert.ToInt32(Math.Pow(ParentalGenotypes.Count, OffspringsIds.Count)));
            GenerateOffspringGenotypes(OffspringsGenotypes, ParentalGenotypes, OffspringsIds.Count, new List<Genotype>());
            if (OffspringsGenotypes.Count > callerParameters.MaxNumOffspringGenotypes)
            {
                OffspringsGenotypes.Shuffle();
                OffspringsGenotypes = OffspringsGenotypes.Take(callerParameters.MaxNumOffspringGenotypes).ToList();
            }
            TransitionMatrix = GetTransitionMatrix(callerParameters.MaximumCopyNumber);
        }


        public List<Genotype> GenerateParentalGenotypes(int numCnStates)
        {
            var genotypes = new List<Genotype>();
            for (int cn = 0; cn < numCnStates; cn++)
            {
                for (int gt = 0; gt <= cn; gt++)
                {
                    genotypes.Add(new Genotype(gt, cn - gt));
                }
            }
            return genotypes;
        }

        public void GenerateOffspringGenotypes(List<List<Genotype>> offspringGenotypes, List<Genotype> genotypeSet,
            int nOffsprings, List<Genotype> partialGenotypes)
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

        public double[][] GetTransitionMatrix(int numCnStates)
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
