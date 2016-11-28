using System;
using System.CodeDom;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using CanvasCommon;
using MathNet.Numerics.Distributions;
using System.Threading.Tasks;
using Isas.SequencingFiles;

namespace CanvasPedigreeCaller
{
    class CanvasPedigreeCaller
    {
        #region Members
        // Static:
        static int MaximumCopyNumber = 5;
        private const int DefaultAlleleDensityThreshold = 1000;
        private const int DefaultAlleleCountThreshold = 4;
        private const double DeNovoRate = 0.00001;

        public int QualityFilterThreshold { get; set; } = 10;
        #endregion

        internal int CallVariants(List<string> variantFrequencyFiles, List<string> segmentFiles, List<string> outVcfFiles, string ploidyBedPath, string referenceFolder, List<string> sampleNames, string pedigreeFile)
        {
            // load files
            // initialize data structures and classes
            int fileCounter = 0;
            Dictionary<string, PedigreeMember.Kinship> kinships = ReadPedigreeFile(pedigreeFile);
            List<PedigreeMember> pedigreeMembers = new List<PedigreeMember>();
            foreach (string sampleName in sampleNames)
            {
                PedigreeMember pedigreeMember = new PedigreeMember();
                pedigreeMember.Name = sampleName;
                pedigreeMember.Segments = CanvasSegment.ReadSegments(segmentFiles[fileCounter]);
                pedigreeMember.MeanMafCoverage = CanvasIO.LoadFrequencies(variantFrequencyFiles[fileCounter], pedigreeMember.Segments);
                foreach (CanvasSegment segment in pedigreeMember.Segments)
                    if (segment.Alleles.Counts.Count > DefaultAlleleCountThreshold)
                        segment.Alleles.SetMedianCounts();
                pedigreeMember.Variance = Math.Pow(Utilities.StandardDeviation(pedigreeMember.Segments.Select(x => x.MedianCount).ToArray()), 2);
                pedigreeMember.MafVariance = Math.Pow(Utilities.StandardDeviation(pedigreeMember.Segments.Where(x => x.Alleles.TotalCoverage.Count > 0).Select(x => x.Alleles.TotalCoverage.Average()).ToArray()), 2);
                pedigreeMember.MeanCoverage = pedigreeMember.Segments.Select(x => x.MedianCount).Average();
                pedigreeMember.MaxCoverage = Convert.ToInt32(pedigreeMember.Segments.Select(x => x.MedianCount).Max() + 10);
                pedigreeMember.Ploidy = PloidyInfo.LoadPloidyFromVcfFile(ploidyBedPath, pedigreeMember.Name);
                pedigreeMember.Kin = kinships[pedigreeMember.Name] == PedigreeMember.Kinship.Parent ?
                    PedigreeMember.Kinship.Parent : PedigreeMember.Kinship.Offspring;
                if (kinships[pedigreeMember.Name] == PedigreeMember.Kinship.Proband)
                    pedigreeMember.Kin = PedigreeMember.Kinship.Proband;
                pedigreeMembers.Add(pedigreeMember);
                fileCounter++;
            }

            var numberOfSegments = pedigreeMembers.First().Segments.Count;
            List<GenomicInterval> segmentIntervals = GetParallelIntevals(numberOfSegments, 1);
            
            List<PedigreeMember> parents = GetParents(pedigreeMembers);
            List<PedigreeMember> offsprings = GetChildren(pedigreeMembers);
            double[][] transitionMatrix = GetTransitionMatrix(MaximumCopyNumber);
            List<Tuple<int, int>> parentalGenotypes = GenerateParentalGenotypes(MaximumCopyNumber);
            List<List<Tuple<int, int>>> offspringsGenotypes = new List<List<Tuple<int, int>>>(Convert.ToInt32(Math.Pow(parentalGenotypes.Count, offsprings.Count)));
            GenerateOffspringGenotypes(offspringsGenotypes, parentalGenotypes, offsprings.Count, new List<Tuple<int, int>>());
            
            foreach (PedigreeMember parent in parents)
                parent.CnModel = new CopyNumberModel(MaximumCopyNumber, parent.MeanCoverage/2.0, parent.MeanMafCoverage/2.0, parent.Variance, parent.MafVariance, parent.MaxCoverage);
            foreach (PedigreeMember offspring in offsprings)
                offspring.CnModel = new CopyNumberModel(MaximumCopyNumber, offspring.MeanCoverage/2.0, offspring.MeanMafCoverage/2.0, offspring.Variance, offspring.MafVariance, offspring.MaxCoverage);

            Parallel.ForEach(
                segmentIntervals,
                new ParallelOptions
                {
                    MaxDegreeOfParallelism = 1,
                    TaskScheduler = TaskScheduler.Default
                },
                interval =>
                {
                    Console.WriteLine($"{DateTime.Now} Launching SPW task for segment {interval.Start} - {interval.End}");
                    var segmentIndex = 0;
                    while (segmentIndex < numberOfSegments)
                    {
                        if (segmentIndex > interval.Start && segmentIndex < interval.End)
                        {
                            var alleleCounts    = pedigreeMembers.Select(x => x.Segments[segmentIndex].Alleles.Counts.Count);
                            var enumerable = alleleCounts as int[] ?? alleleCounts.ToArray();
                            var alleleDensity   = enumerable.Average()/pedigreeMembers.First().Segments[segmentIndex].Length;
                            var useCnLikelihood = enumerable.Select(x => x > DefaultAlleleCountThreshold).Any(c => c == false) && 
                                alleleDensity < DefaultAlleleDensityThreshold;
                            CopyNumberDistribution copyNumberLikelihoods = useCnLikelihood ? 
                            MaximalCnLikelihood(parents, offsprings, segmentIndex, transitionMatrix, offspringsGenotypes) : 
                            MaximalGtLikelihood(parents, offsprings, segmentIndex, parentalGenotypes, offspringsGenotypes);
                            EstimateQScores(parents, offsprings, segmentIndex, copyNumberLikelihoods);
                        }

                        segmentIndex++;
                    }
                    Console.WriteLine($"{DateTime.Now} Finished SPW task for segment {interval.Start} - {interval.End}");
                });

            int pedigreeMemberIndex = 0;
            foreach (var pedigreeMember in pedigreeMembers)
            {
                CanvasSegment.MergeSegments(ref pedigreeMember.Segments);
                CanvasSegment.WriteSegments(outVcfFiles[pedigreeMemberIndex], pedigreeMember.Segments,
                    pedigreeMember.MeanCoverage, referenceFolder, pedigreeMember.Name, null, null, QualityFilterThreshold);
                pedigreeMemberIndex++;
            }
            return 0;
        }

        private void EstimateQScores(List<PedigreeMember> parents, List<PedigreeMember> offsprings, int segmentIndex,
            CopyNumberDistribution copyNumberLikelihoods)
        {
            var cnStates = GetCnStates(parents, offsprings, segmentIndex);
            var names = parents.Concat(offsprings).Select(x => x.Name).ToList();
            var probands = GetProbands(offsprings);
            var singleSampleQualityScores = GetSingleSampleQualityScores(copyNumberLikelihoods, cnStates, names);
            if (parents.First().Segments[segmentIndex].Chr == "chr1" &&
                parents.First().Segments[segmentIndex].Begin > 229812690 &&
                parents.First().Segments[segmentIndex].Begin < 229812699)
            {
                var cnMarginalProbabilities = copyNumberLikelihoods.GetMarginalProbability(cnStates.Count, MaximumCopyNumber, names[0]);
                foreach (var prob in cnMarginalProbabilities)
                {
                    Console.WriteLine($"prob = {prob}");
                }
            }
            var probandIndex = names.IndexOf(probands.First().Name);
            var parent1Index = names.IndexOf(parents.First().Name);
            var parent2Index = names.IndexOf(parents.Last().Name);

            if (cnStates[probandIndex] != 2 && cnStates[parent1Index] == 2 && cnStates[parent2Index] == 2 && 
                singleSampleQualityScores[probandIndex] > 20)
            {
                var deNovoQualityScore = GetDeNovoQualityScore(parents, copyNumberLikelihoods, probands.First().Name,
                    cnStates[probandIndex], singleSampleQualityScores[probandIndex]);
                probands.First().Segments[segmentIndex].DQScore = deNovoQualityScore;
            }
            var counter = 0;
            foreach (PedigreeMember sample in parents.Concat(offsprings))
            {
                sample.Segments[segmentIndex].QScore = singleSampleQualityScores[counter];
                counter++;
            }
        }

        private static List<int> GetCnStates(IEnumerable<PedigreeMember> parents, IEnumerable<PedigreeMember> offsprings, int segmentIndex)
        {
            return parents.Concat(offsprings).Select(x=>Math.Min(x.Segments[segmentIndex].CopyNumber, MaximumCopyNumber-1)).ToList();
        }

        private static List<PedigreeMember> GetChildren(IEnumerable<PedigreeMember> pedigreeMembers)
        {
            return pedigreeMembers.Select(x => x).Where(x => x.Kin != PedigreeMember.Kinship.Parent).ToList();
        }

        private static List<PedigreeMember> GetParents(IEnumerable<PedigreeMember> pedigreeMembers)
        {
            return pedigreeMembers.Select(x => x).Where(x => x.Kin == PedigreeMember.Kinship.Parent).ToList();
        }

        private static List<PedigreeMember> GetProbands(IEnumerable<PedigreeMember> pedigreeMembers)
        {
            return pedigreeMembers.Select(x => x).Where(x => x.Kin == PedigreeMember.Kinship.Proband).ToList();
        }

        public static int AggregateVariantCoverage(ref List<CanvasSegment> segments)
        {
            var variantCoverage = segments.SelectMany(segment => segment.Alleles.TotalCoverage).ToList();
            return variantCoverage.Any() ? Utilities.Median(variantCoverage) : 0;
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

        /// <summary>
        /// Calculates maximal likelihood for segments without SNV allele ratios. Updated CanvasSegment CopyNumber only. 
        /// </summary>
        /// <param name="parents"></param>
        /// <param name="children"></param>
        /// <param name="segmentPosition"></param>
        /// <param name="transitionMatrix"></param>
        public CopyNumberDistribution MaximalCnLikelihood(List<PedigreeMember> parents, List<PedigreeMember> children, int segmentPosition, double[][] transitionMatrix, List<List<Tuple<int, int>>> offspringsGenotypes)
        {
            double maximalLikelihood;
            double marginals;
            InitializeLikelihood(out maximalLikelihood, out marginals, segmentPosition, parents, children);
            var parent1Likelihood = parents.First().CnModel.GetCnLikelihood(Math.Min(parents.First().GetCoverage(segmentPosition), parents.First().MeanCoverage*3.0));
            var parent2Likelihood = parents.Last().CnModel.GetCnLikelihood(Math.Min(parents.Last().GetCoverage(segmentPosition), parents.Last().MeanCoverage*3.0));

            if (parent1Likelihood.Count != parent2Likelihood.Count)
                throw new ArgumentException("Both parents should have the same number of CN states");
            int nCopies = MaximumCopyNumber;
            List<string> names = parents.Select(x => x.Name).Union(children.Select(x => x.Name)).ToList();
            CopyNumberDistribution density = new CopyNumberDistribution(nCopies, names);

            for (int cn1 = 0; cn1 < nCopies; cn1++)
            {
                for (int cn2 = 0; cn2 < nCopies; cn2++)
                {
                    foreach (var offspringGtStates in offspringsGenotypes)
                    {
                        double currentLikelihood = parent1Likelihood[cn1] * parent2Likelihood[cn2];
                        int counter = 0;
                        foreach (PedigreeMember child in children)
                        {
                            var modelIndex = Math.Min(offspringGtStates[counter].Item1 + offspringGtStates[counter].Item2, MaximumCopyNumber - 1);
                            currentLikelihood *= transitionMatrix[cn1][offspringGtStates[counter].Item1] * 
                                transitionMatrix[cn2][offspringGtStates[counter].Item2] *
                                child.CnModel.GetCnLikelihood(child.GetCoverage(segmentPosition))[modelIndex];
                            counter++;
                            
                        }
                        int[] copyNumberIndices = { cn1, cn2 };
                        var index = copyNumberIndices.Concat(offspringGtStates.Select(x => x.Item1 + x.Item2)).ToArray();
                        density.SetJointProbability(Math.Max(currentLikelihood, density.GetJointProbability(index)), index);

                        currentLikelihood = Double.IsNaN(currentLikelihood) || Double.IsInfinity(currentLikelihood) ? 0 : currentLikelihood;

                        if (currentLikelihood > maximalLikelihood)
                        {
                            maximalLikelihood = currentLikelihood;
                            parents.First().Segments[segmentPosition].CopyNumber = cn1;
                            parents.Last().Segments[segmentPosition].CopyNumber  = cn2;
                            counter = 0;
                            foreach (PedigreeMember child in children)
                            {
                                child.Segments[segmentPosition].CopyNumber = offspringGtStates[counter].Item1 + offspringGtStates[counter].Item2;
                                counter++;
                            }
                        }
                        marginals += maximalLikelihood;
                    }
                }
            }     
            return density;
        }

        /// <summary>
        /// Calculates maximal likelihood for segments with SNV allele  counts. Updated CanvasSegment CopyNumber and MajorChromosomeCount.
        /// </summary>
        /// <param name="parents"></param>
        /// <param name="children"></param>
        /// <param name="segmentPosition"></param>
        public CopyNumberDistribution MaximalGtLikelihood(List<PedigreeMember> parents, List<PedigreeMember> children, int segmentPosition, List<Tuple<int, int>> parentalGenotypes, List<List<Tuple<int, int>>> offspringsGenotypes)
        {
            double maximalLikelihood;
            double marginals;
            InitializeLikelihood(out maximalLikelihood, out marginals, segmentPosition, parents, children);
            var parent1Likelihood = parents.First().CnModel.GetGtLikelihood(parents.First().GetAlleleCounts(segmentPosition));
            var parent2Likelihood = parents.Last().CnModel.GetGtLikelihood(parents.Last().GetAlleleCounts(segmentPosition));
            int nCopies = MaximumCopyNumber;
            List<string> names = parents.Select(x => x.Name).Union(children.Select(x => x.Name)).ToList();
            CopyNumberDistribution density = new CopyNumberDistribution(nCopies, names);

            foreach (var parent1GtStates in parentalGenotypes)
            {
                foreach (var parent2GtStates in parentalGenotypes)
                {
                    foreach (var offspringGtStates in offspringsGenotypes)
                    {
                        var currentLikelihood = parent1Likelihood[parent1GtStates.Item1][parent1GtStates.Item2] *
                        parent2Likelihood[parent2GtStates.Item1][parent2GtStates.Item2];
                        int counter = 0;
                        foreach (PedigreeMember child in children)
                        {
                            currentLikelihood *= GetTransition(parent1GtStates.Item1, parent1GtStates.Item2, offspringGtStates[counter].Item1, offspringGtStates[counter].Item2) *
                                                 GetTransition(parent2GtStates.Item1, parent2GtStates.Item2, offspringGtStates[counter].Item1, offspringGtStates[counter].Item2) *
                                                 child.CnModel.GetGtLikelihood(child.GetAlleleCounts(segmentPosition))[offspringGtStates[counter].Item1][offspringGtStates[counter].Item2];
                            counter++;
                        }

                        int[] copyNumberIndices = { parent1GtStates.Item1 + parent1GtStates.Item2, parent2GtStates.Item1 + parent2GtStates.Item2 };
                        var index = copyNumberIndices.Concat(offspringGtStates.Select(x => x.Item1 + x.Item2)).ToArray();
                        if (density.GetJointProbability(index) < currentLikelihood)
                            density.SetJointProbability(Math.Max(currentLikelihood, Double.MinValue), index);

                        currentLikelihood = Double.IsNaN(currentLikelihood) || Double.IsInfinity(currentLikelihood) ? 0 : currentLikelihood;

                        if (currentLikelihood > maximalLikelihood)
                        {
                            maximalLikelihood = currentLikelihood;
                            parents.First().Segments[segmentPosition].CopyNumber = parent1GtStates.Item1 + parent1GtStates.Item2;
                            parents.First().Segments[segmentPosition].MajorChromosomeCount = Math.Max(parent1GtStates.Item1, parent1GtStates.Item2);

                            parents.Last().Segments[segmentPosition].CopyNumber = parent2GtStates.Item1 + parent2GtStates.Item2;
                            parents.Last().Segments[segmentPosition].MajorChromosomeCount = Math.Max(parent2GtStates.Item1, parent2GtStates.Item2);

                            counter = 0;
                            foreach (PedigreeMember child in children)
                            {
                                child.Segments[segmentPosition].CopyNumber = parent2GtStates.Item1 + offspringGtStates[counter].Item2;
                                child.Segments[segmentPosition].MajorChromosomeCount = Math.Max(offspringGtStates[counter].Item1, offspringGtStates[counter].Item2);
                                counter++;
                            }
                        }
                        marginals += maximalLikelihood;
                    }
                }
            }
            return density;
        }

        private static void InitializeLikelihood(out double maximalLikelihood, out double marginals, int segmentPosition, List<PedigreeMember> parents, List<PedigreeMember> children)
        {
            maximalLikelihood = 0;
            marginals = 0;
            int defaultCn = 2;

            parents.First().Segments[segmentPosition].CopyNumber = defaultCn;
            parents.Last().Segments[segmentPosition].CopyNumber = defaultCn;
            foreach (PedigreeMember child in children)
                child.Segments[segmentPosition].CopyNumber = defaultCn;
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

        public void GenerateOffspringGenotypes(List<List<Tuple<int, int>>> offspringGenotypes, List<Tuple<int, int>> genotypeSet, int nOffsprings, List<Tuple<int, int>> partialGenotypes)
        {      

            if (nOffsprings > 0)
            {
                foreach (Tuple<int, int> genotype in genotypeSet)
                {
                    GenerateOffspringGenotypes(offspringGenotypes, genotypeSet, nOffsprings - 1, partialGenotypes.Concat(new List<Tuple<int, int>> { genotype }).ToList());
                }
            }
            if (nOffsprings == 0)
            {
                offspringGenotypes.Add(partialGenotypes);
            }                  
        }

        private List<GenomicInterval> GetParallelIntevals(int nSegments, int nCores)
        {
            List<GenomicInterval> intevals = new List<GenomicInterval>();
            int step = nSegments/nCores;
            intevals.Add(new GenomicInterval(0, step));
            int cumSum = step + 1;
            while (cumSum + step + 1 < nSegments-1)
            {
                intevals.Add(new GenomicInterval(cumSum, cumSum + step));
                cumSum += step+1;
            }
            intevals.Add(new GenomicInterval(cumSum, nSegments-1));
            return intevals;
        }

        public double GetTransition(int gt1Parent, int gt2Parent, int gt1Offspring, int gt2Offspring)
        {
            if (gt1Parent == gt1Offspring || gt1Parent == gt2Offspring ||
                gt2Parent == gt1Offspring || gt2Parent == gt2Offspring)
                return 0.5;
            return DeNovoRate;
        }


        public Dictionary<string, PedigreeMember.Kinship> ReadPedigreeFile(string pedigreeFile)
        {
            Dictionary<string, PedigreeMember.Kinship> kinships = new Dictionary<string, PedigreeMember.Kinship>();
            using (StreamReader reader = new StreamReader(pedigreeFile))
            {
                string row;
                while ((row = reader.ReadLine()) != null)
                {
                    string[] fields = row.Split('\t');
                    string maternalId = fields[2];
                    string paternallId = fields[3];
                    string proband = fields[5];
                    if (maternalId == "0" && paternallId == "0")
                        kinships.Add(fields[1], PedigreeMember.Kinship.Parent);
                    else if (proband == "affected")
                        kinships.Add(fields[1], PedigreeMember.Kinship.Proband);
                    else
                        kinships.Add(fields[1], PedigreeMember.Kinship.Offspring);
                }
            }
            return kinships;
        }

        public List<double> GetSingleSampleQualityScores(CopyNumberDistribution density, List<int> cnStates, List<string> sampleNames)
        {
            var singleSampleQualityScores = new List<double>();
            if (density.Count != cnStates.Count)
                throw new ArgumentException("Size of CopyNumberDistribution should be equal to number of CN states");
            for (int index = 0; index < sampleNames.Count; index++)
            {
                string sampleName = sampleNames[index];
                var cnMarginalProbabilities = density.GetMarginalProbability(cnStates.Count, MaximumCopyNumber, sampleName);
                double normalizationConstant = cnMarginalProbabilities.Sum();
                var qscore = -10.0*Math.Log10((normalizationConstant - cnMarginalProbabilities[cnStates[index]]) / normalizationConstant);
                singleSampleQualityScores.Add(qscore);
            }
            return singleSampleQualityScores;
        }

        public double GetDeNovoQualityScore(List<PedigreeMember> parents, CopyNumberDistribution density, string sampleName, int sampleValue, double sampleProbability)
        {
            int nSamples = density.Count;
            var cnMarginalProbabilities = density.GetMarginalProbability(nSamples, MaximumCopyNumber, sampleName);
            density.SetConditionalProbability(density.Count, MaximumCopyNumber, sampleName,  sampleValue, cnMarginalProbabilities[sampleValue]);
            var parentNames = parents.Select(x => x.Name).ToList();
            var parentCopyNumber = 2;
            var firstParentMarginalProbabilities = density.GetMarginalProbability(nSamples, MaximumCopyNumber, parentNames.First());
            var secondParentMarginalProbabilities = density.GetMarginalProbability(nSamples, MaximumCopyNumber, parentNames.Last());
            var qscorefirstParent = -10.0 * Math.Log10((firstParentMarginalProbabilities.Sum() - firstParentMarginalProbabilities[parentCopyNumber]) / 
                firstParentMarginalProbabilities.Sum());
            var qscoresecondParent = -10.0 * Math.Log10((secondParentMarginalProbabilities.Sum() - secondParentMarginalProbabilities[parentCopyNumber]) /
    secondParentMarginalProbabilities.Sum());
            return (qscorefirstParent + qscoresecondParent) / 2.0;
        }
}
}
