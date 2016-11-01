using System;
using System.Collections.Generic;
using System.Diagnostics;
using Illumina.Common;
using System.IO;
using System.Linq;
using CanvasCommon;
using MathNet.Numerics.Distributions;
using System.Threading.Tasks;

namespace CanvasPedigreeCaller
{
    class CanvasPedigreeCaller
    {
        #region Members
        // Static:
        static int MaximumCopyNumber = 5;
        public int QualityFilterThreshold { get; set; } = 10;
        #endregion

        internal int CallVariants(List<string> variantFrequencyFiles, List<string> segmentFiles, List<string> outDirs, List<string> ploidyBedPaths, string referenceFolder, List<string> sampleNames, string pedigreeFile)
        {
            // load files
            // initialize data structures and classes
            int fileCounter = 0;
            Dictionary<string, Kinship> kinships = ReadPedigreeFile(pedigreeFile);
            List<PedigreeMember> pedigreeMembers = new List<PedigreeMember>();
            foreach (string sampleName in sampleNames)
            {
                PedigreeMember pedigreeMember = new PedigreeMember();
                pedigreeMember.Name = sampleName;
                pedigreeMember.Segments = CanvasSegment.ReadSegments(segmentFiles[fileCounter]);
                pedigreeMember.MeanMafCoverage = CanvasIO.LoadVariantFrequencies(variantFrequencyFiles[fileCounter], pedigreeMember.Segments);
                pedigreeMember.MeanCoverage = pedigreeMember.Segments.Select(x => x.Counts.Average()).Average();
                pedigreeMember.Ploidy = PloidyInfo.LoadPloidyFromBedFile(ploidyBedPaths[fileCounter]);
                pedigreeMembers.Add(pedigreeMember);
            }
            Pedigree<PedigreeMember> pedigree = new Pedigree<PedigreeMember>(pedigreeMembers[0], true);
            foreach (PedigreeMember pedigreeMember in pedigreeMembers)
            {
                if (kinships[pedigreeMember.Name] == Kinship.Offspring)
                    pedigree.AddChildren(pedigreeMember);
                else
                    pedigree.AddParent(pedigreeMember);
            }

            var numberOfSegments = pedigreeMembers.First().Segments.Count;
            List<string> chrmosomes = new List<string>();
            foreach (CanvasSegment segment in pedigreeMembers.First().Segments)
                if (!chrmosomes.Contains(segment.Chr))
                    chrmosomes.Add(segment.Chr);

            List<PedigreeMember> parents = pedigree.GetParents();
            List<PedigreeMember> offsprings = pedigree.GetChildren();
            double[][] transitionMatrix = GetTransitionMatrix(MaximumCopyNumber);
            List<Tuple<int, int>> parentalGenotypes = GenerateParentalGenotypes(MaximumCopyNumber);
            List<List<Tuple<int, int>>> offspringsGenotypes = GenerateOffspringGenotypes(parentalGenotypes, offsprings.Count);
            

            foreach (PedigreeMember parent in parents)
                parent.CnModel = new CopyNumberModel(MaximumCopyNumber, parent.MeanCoverage/2.0, parent.MeanMafCoverage/2.0, parent.Variance);
            foreach (PedigreeMember offspring in offsprings)
                offspring.CnModel = new CopyNumberModel(MaximumCopyNumber, offspring.MeanCoverage/2.0, offspring.MeanMafCoverage/2.0, offspring.Variance);

            Parallel.ForEach(
                chrmosomes,
                new ParallelOptions
                {
                    MaxDegreeOfParallelism = Environment.ProcessorCount,
                    TaskScheduler = TaskScheduler.Default
                },
                chr =>
                {
                    Console.WriteLine($"{DateTime.Now} Launching SPW task for chromosome {chr}");
                    var segmentIndex = 0;
                    while (segmentIndex < numberOfSegments)
                    {
                        if (pedigreeMembers.First().Segments[segmentIndex].Chr == chr)
                        {
                            if (pedigreeMembers.Select(x => x.Segments[segmentIndex].VariantAlleleCounts.Count > 0).Any(c => c == false))
                                MaximalCnLikelyhood(parents, offsprings, segmentIndex, transitionMatrix, offspringsGenotypes);
                            else
                                MaximalGtLikelyhood(parents, offsprings, segmentIndex, parentalGenotypes, offspringsGenotypes);
                        }
                        segmentIndex++;
                    }
                    Console.WriteLine($"{DateTime.Now} Finished SPW task for chromosome {chr}");
                });


            for (int pedigreeMemberIndex = 0; pedigreeMemberIndex < pedigreeMembers.Count; pedigreeMemberIndex++)
            {
                var outFile = Path.Combine(outDirs[pedigreeMemberIndex], pedigreeMembers[pedigreeMemberIndex].Name + "CNV.vcf");
                CanvasSegment.MergeSegments(ref pedigreeMembers[pedigreeMemberIndex].Segments);
                CanvasSegment.WriteCoveragePlotData(pedigreeMembers[pedigreeMemberIndex].Segments,
                    pedigreeMembers[pedigreeMemberIndex].MeanCoverage * 2.0, pedigreeMembers[pedigreeMemberIndex].Ploidy, outFile, referenceFolder);
            }
            return 0;
        }


        public static int AggregateVariantCoverage(ref List<CanvasSegment> segments)
        {
            var variantCoverage = segments.SelectMany(segment => segment.VariantTotalCoverage).ToList();
            return variantCoverage.Any() ? Utilities.Median(variantCoverage) : 0;
        }


        double[][] GetTransitionMatrix(int numCnStates)
        {
            double[][] transitionMatrix = CanvasCommon.Utilities.MatrixCreate(numCnStates, numCnStates);
            for (int cn = 0; cn < numCnStates; cn++)
            {
                var gtLikelyhood = new Poisson(Math.Max(cn / 2.0, 0.1));
                for (int gt = 0; gt < cn; gt++)
                    transitionMatrix[cn][gt] = gtLikelyhood.Probability(gt);
            }
            return transitionMatrix;
        }

        /// <summary>
        /// Calculates maximal likelihood for segments without SNV allele ratios. Updated CanvasSegment CopyNumber only. 
        /// </summary>
        /// <param name="parents"></param>
        /// <param name="children"></param>
        /// <param name="segmentPosition"></param>
        /// <param name="transitionMatrix"></param>
        public void MaximalCnLikelyhood(List<PedigreeMember> parents, List<PedigreeMember> children, int segmentPosition, double[][] transitionMatrix, List<List<Tuple<int, int>>> offspringsGenotypes)
        {
            double maximalLikelyhood = Double.MinValue;
            double marginals = 0;
            List<int> parentCnStates = new List<int>();
            var parent1Likelihood = parents[0].CnModel.GetCnLikelihood(parents[0].GetCoverage(segmentPosition));
            var parent2Likelihood = parents[1].CnModel.GetCnLikelihood(parents[1].GetCoverage(segmentPosition));
            for (int cn1 = 0; cn1 < parent1Likelihood.Count; cn1++)
            {
                parentCnStates.Add(cn1);
                for (int cn2 = 0; cn2 < parent2Likelihood.Count; cn2++)
                {
                    parentCnStates.Add(cn2);
                    for (int gt1 = 0; gt1 < cn1; gt1++)
                    {
                        for (int gt2 = 0; gt2 < cn2; gt2++)
                        {
                            foreach (List<Tuple<int, int>> offspringGtStates in offspringsGenotypes)
                            {
                                double currentLikelihood = parent1Likelihood[cn1]*parent2Likelihood[cn2];
                                int counter = 0;
                                foreach (PedigreeMember child in children)
                                {
                                    currentLikelihood *= transitionMatrix[cn1][offspringGtStates[counter].Item1] *transitionMatrix[cn2][offspringGtStates[counter].Item2] *
                                                         child.CnModel.GetCnLikelihood(child.GetCoverage(segmentPosition))
                                                             [Math.Min(offspringGtStates[counter].Item1 + offspringGtStates[counter].Item2, MaximumCopyNumber - 1)];
                                    counter++;
                                }
                                if (currentLikelihood > maximalLikelyhood)
                                {
                                    maximalLikelyhood = currentLikelihood;
                                     counter = 0;
                                    foreach (PedigreeMember parent in parents)
                                    {
                                        parent.Segments[segmentPosition].CopyNumber = parentCnStates[counter];
                                        counter++;
                                    }
                                    counter = 0;
                                    foreach (PedigreeMember child in children)
                                    {
                                        child.Segments[segmentPosition].CopyNumber = offspringGtStates[counter].Item1 + offspringGtStates[counter].Item2;
                                        counter++;
                                    }
                                }
                                marginals += maximalLikelyhood;
                            }
                        }
                    }
                    parentCnStates.RemoveAt(parentCnStates.Count - 1);
                }
                parentCnStates.Clear();
            }

            foreach (PedigreeMember parent in parents)
                parent.Segments[segmentPosition].QScore = marginals == 0 ? 0 : maximalLikelyhood / marginals;

            foreach (PedigreeMember offring in children)
                offring.Segments[segmentPosition].QScore = marginals == 0 ? 0 : maximalLikelyhood / marginals;
        }

        /// <summary>
        /// Calculates maximal likelihood for segments with SNV allele  counts. Updated CanvasSegment CopyNumber and MajorChromosomeCount.
        /// </summary>
        /// <param name="parents"></param>
        /// <param name="children"></param>
        /// <param name="segmentPosition"></param>
        public void MaximalGtLikelyhood(List<PedigreeMember> parents, List<PedigreeMember> children, int segmentPosition, List<Tuple<int, int>> parentalGenotypes, List<List<Tuple<int, int>>> offspringsGenotypes)
        {
            double maximalLikelyhood = Double.MinValue;
            double marginals = 0;
            var parent1Likelihood = parents[0].CnModel.GetGtLikelihood(parents[0].GetAlleleCounts(segmentPosition));
            var parent2Likelihood = parents[1].CnModel.GetGtLikelihood(parents[1].GetAlleleCounts(segmentPosition));
            int maxCns = parent1Likelihood.Length;
            List<Tuple<int, int>> parentGtStates = new List<Tuple<int, int>>();
            foreach (Tuple<int, int> parent1GtStates in parentalGenotypes)
            {
                parentGtStates.Add(parent1GtStates);
                foreach (Tuple<int, int> parent2GtStates in parentalGenotypes)
                {
                    parentGtStates.Add(parent2GtStates);
                    foreach (List<Tuple<int, int>> offspringGtStates in offspringsGenotypes)
                    {
                        double currentLikelihood = parent1Likelihood[parent1GtStates.Item1][parent1GtStates.Item2] *
                        parent2Likelihood[parent2GtStates.Item1][parent2GtStates.Item2];
                        int counter = 0;
                        foreach (PedigreeMember child in children)
                        {
                            currentLikelihood *= getTransition(parent1GtStates.Item1, parent1GtStates.Item2, offspringGtStates[counter].Item1, offspringGtStates[counter].Item2) *
                                                 getTransition(parent2GtStates.Item1, parent2GtStates.Item2, offspringGtStates[counter].Item1, offspringGtStates[counter].Item2) *
                                                 child.CnModel.GetGtLikelihood(child.GetAlleleCounts(segmentPosition))[offspringGtStates[counter].Item1][offspringGtStates[counter].Item2];
                            counter++;
                        }
                        if (currentLikelihood > maximalLikelyhood)
                        {
                            maximalLikelyhood = currentLikelihood;
                            counter = 0;
                            foreach (PedigreeMember parent in parents)
                            {
                                parent.Segments[segmentPosition].CopyNumber = parentGtStates[counter].Item1 +
                                                                              parentGtStates[counter].Item2;
                                parent.Segments[segmentPosition].MajorChromosomeCount =
                                    Math.Max(parentGtStates[counter].Item1, parentGtStates[counter].Item2);
                                counter++;
                            }
                            counter = 0;
                            foreach (PedigreeMember child in children)
                            {
                                child.Segments[segmentPosition].CopyNumber = offspringGtStates[counter].Item1 +
                                                                              offspringGtStates[counter].Item2;
                                child.Segments[segmentPosition].MajorChromosomeCount =
                                    Math.Max(offspringGtStates[counter].Item1, offspringGtStates[counter].Item2);
                                counter++;
                            }
                        }
                        marginals += maximalLikelyhood;
                    }
                    parentGtStates.RemoveAt(parentGtStates.Count - 1);
                }
                parentGtStates.Clear();
            }

            foreach (PedigreeMember parent in parents)
                parent.Segments[segmentPosition].QScore = marginals == 0 ? 0 : maximalLikelyhood / marginals;

            foreach (PedigreeMember offring in children)
                offring.Segments[segmentPosition].QScore = marginals == 0 ? 0 : maximalLikelyhood / marginals;
        }


        public List<Tuple<int, int>> GenerateParentalGenotypes(int numCnStates)
        {
            List<Tuple<int, int>> genotypes = new List<Tuple<int, int>>();
            for (int cn = 0; cn < numCnStates; cn++)
            {
                for (int gt = 0; gt <= cn; gt++)
                {
                    genotypes.Add(new Tuple<int, int>(gt, cn - gt));
                }
            }
            return genotypes;
        }

        ///
        public List<List<Tuple<int, int>>> GenerateOffspringGenotypes(List<Tuple<int, int>> genotypes, int nOffsprings)
        {
            List<List<Tuple<int, int>>> offspringsGenotypes = new List<List<Tuple<int, int>>>();
            int nOffspringCounter = nOffsprings;
            int iterationCounter = nOffsprings;

            while (nOffspringCounter > 0)
            {
                foreach (Tuple<int, int> genotype in genotypes)
                {
                    offspringsGenotypes.Add(new List<Tuple<int, int>> {genotype});
                }
                nOffspringCounter--;
            }
            iterationCounter--;

            while (iterationCounter > 0)
            {
                nOffspringCounter = nOffsprings;
                int genotypeCounter = 0;
                while (nOffspringCounter > 0)
                {
                    foreach (Tuple<int, int> genotype in genotypes)
                    {
                        offspringsGenotypes[genotypeCounter].Add(genotype);
                        genotypeCounter++;
                    }
                    nOffspringCounter--;
                }
                iterationCounter--;
            }
            return offspringsGenotypes;
        }


        public double getTransition(int gt1Parent, int gt2Parent, int gt1Offspring, int gt2Offspring)
        {
            if (gt1Parent == gt1Offspring || gt1Parent == gt2Offspring ||
                gt2Parent == gt1Offspring || gt2Parent == gt2Offspring)
                return 0.5;
            return 0.00001;
        }

        public enum Kinship
        {
            Parent, Offspring
        }
        public Dictionary<string, Kinship> ReadPedigreeFile(string pedigreeFile)
        {
            Dictionary<string, Kinship> kinships = new Dictionary<string, Kinship>();
            using (StreamReader reader = new StreamReader(pedigreeFile))
            {
                string row;
                while ((row = reader.ReadLine()) != null)
                {
                    string[] fields = row.Split('\t');
                    string maternalId = fields[2];
                    string paternallId = fields[3];
                    if (maternalId == "0" && paternallId == "0")
                        kinships.Add(fields[1], Kinship.Parent);
                    else
                        kinships.Add(fields[1], Kinship.Offspring);
                }
            }
            return kinships;
        }
    }
}
