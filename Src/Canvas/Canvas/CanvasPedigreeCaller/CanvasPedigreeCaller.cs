using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.Linq;
using System.Linq.Expressions;
using System.Xml.Schema;
using CanvasCommon;
using Illumina.Common;
using MathNet.Numerics.Distributions;
using System.Threading.Tasks;

namespace CanvasPedigreeCaller
{
    class CanvasPedigreeCaller
    {
        #region Members
        // Static:
        static private int MaximumCopyNumber = 10;
        public int QualityFilterThreshold { get; set; } = 10;
        #endregion


        internal int CallVariants(List<string> variantFrequencyFiles, List<string> segmentFiles, string outDir, List<string> ploidyBedPaths, string referenceFolder, List<string> sampleNames, string pedigreeFile)
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
                pedigreeMember.MeanCoverage = CanvasIO.LoadVariantFrequencies(variantFrequencyFiles[fileCounter], pedigreeMember.Segments);
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
            double[][] transitionMatrix = GetTransitionMatrix(numberOfSegments);

            foreach (PedigreeMember parent in parents)
                parent.CnModel = new CopyNumberModel(numberOfSegments, parent.MeanCoverage, parent.Variance);
            foreach (PedigreeMember offspring in offsprings)
                offspring.CnModel = new CopyNumberModel(numberOfSegments, offspring.MeanCoverage, offspring.Variance);

            var cts = new System.Threading.Tasks.CancellationTokenSource();
            Parallel.ForEach(
                chrmosomes,
                new ParallelOptions
                {
                    CancellationToken = cts.Token,
                    MaxDegreeOfParallelism = Environment.ProcessorCount,
                    TaskScheduler = TaskScheduler.Default
                },
                chr =>
                {
                    var segmentIndex = 0;                   
                    while (segmentIndex < numberOfSegments)
                    {
                        if (parents.First().Segments[segmentIndex].Chr == chr)
                        {
                            if (pedigreeMembers.Select(x => x.Segments[segmentIndex].VariantAlleleCounts.Empty()).Any(c => c == false))
                                MaximalCnLikelyhood(parents, offsprings, segmentIndex, transitionMatrix);
                            else
                                MaximalGtLikelyhood(parents, offsprings, segmentIndex, transitionMatrix);
                        }
                        segmentIndex++;
                    }
                });


            for (int pedigreeMemberIndex = 0; pedigreeMemberIndex < pedigreeMembers.Count; pedigreeMemberIndex++)
            {
                var outFile = Path.Combine(outDir, pedigreeMembers[pedigreeMemberIndex].Name + "CNV.vcf");
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
                var gtLikelyhood = new Poisson(Math.Round(cn / 2.0));
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
        public void MaximalCnLikelyhood(List<PedigreeMember> parents, List<PedigreeMember> children, int segmentPosition, double[][] transitionMatrix)
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
                            List<int> offspringCnStates = new List<int>();
                            double currentLikelihood = parent1Likelihood[cn1] * parent2Likelihood[cn2];

                            foreach (PedigreeMember child in children)
                            {
                                offspringCnStates.Add(gt1 + gt2);
                                currentLikelihood *= transitionMatrix[cn1][gt1] * transitionMatrix[cn2][gt2] *
                                    child.CnModel.GetCnLikelihood(child.GetCoverage(segmentPosition))[gt1 + gt2];
                            }
                            if (currentLikelihood > maximalLikelyhood)
                            {
                                maximalLikelyhood = currentLikelihood;
                                int counter = 0;
                                foreach (PedigreeMember parent in parents)
                                {
                                    parent.Segments[segmentPosition].CopyNumber = parentCnStates[counter];
                                    counter++;
                                }
                                counter = 0;
                                foreach (PedigreeMember child in children)
                                {
                                    child.Segments[segmentPosition].CopyNumber = offspringCnStates[counter];
                                    counter++;
                                }
                            }
                            marginals += maximalLikelyhood;
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
        /// <param name="transitionMatrix"></param>
        public void MaximalGtLikelyhood(List<PedigreeMember> parents, List<PedigreeMember> children, int segmentPosition, double[][] transitionMatrix)
        {
            double maximalLikelyhood = Double.MinValue;
            double marginals = 0;
            List<Tuple<int,int>> parentGtStates = new List<Tuple<int, int>>();
            var parent1Likelihood = parents[0].CnModel.GetGtLikelihood(parents[0].GetAlleleCounts(segmentPosition));
            var parent2Likelihood = parents[1].CnModel.GetGtLikelihood(parents[1].GetAlleleCounts(segmentPosition));
            int maxCns = parent1Likelihood.Count;

            for (int cn1 = 0; cn1 < maxCns; cn1++)
            {
                parentGtStates.Add(new Tuple<int, int>(cn1, maxCns-cn1));
                for (int cn2 = 0; cn2 < maxCns; cn2++)
                {
                    parentGtStates.Add(new Tuple<int, int>(cn2, maxCns - cn2));
                    for (int gt1 = 0; gt1 < cn1; gt1++)
                    {
                        for (int gt2 = 0; gt2 < cn2; gt2++)
                        {
                            List<Tuple<int, int>> offspringGtStates = new List<Tuple<int, int>>();
                            double currentLikelihood = parent1Likelihood[cn1] * parent2Likelihood[cn2];

                            foreach (PedigreeMember child in children)
                            {
                                offspringGtStates.Add(new Tuple<int, int>(gt1, gt2));
                                currentLikelihood *= transitionMatrix[cn1][gt1] * transitionMatrix[cn2][gt2] *
                                    child.CnModel.GetGtLikelihood(child.GetAlleleCounts(segmentPosition))[gt1 + gt2];
                            }
                            if (currentLikelihood > maximalLikelyhood)
                            {
                                maximalLikelyhood = currentLikelihood;
                                int counter = 0;
                                foreach (PedigreeMember parent in parents)
                                {
                                    parent.Segments[segmentPosition].CopyNumber = parentGtStates[counter].Item1 + parentGtStates[counter].Item1;
                                    parent.Segments[segmentPosition].MajorChromosomeCount = Math.Max(parentGtStates[counter].Item1, parentGtStates[counter].Item1);
                                    counter++;
                                }
                                counter = 0;
                                foreach (PedigreeMember child in children)
                                {
                                    child.Segments[segmentPosition].CopyNumber = Math.Max(offspringGtStates[counter].Item1, offspringGtStates[counter].Item1);
                                    counter++;
                                }
                            }
                            marginals += maximalLikelyhood;
                        }
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
