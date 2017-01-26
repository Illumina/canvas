using System;
using System.CodeDom;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using CanvasCommon;
using MathNet.Numerics.Distributions;
using System.Threading.Tasks;
using System.Xml.Schema;
using Combinatorics.Collections;
using Illumina.Common;
using Isas.SequencingFiles;

namespace CanvasPedigreeCaller
{
    class CanvasPedigreeCaller
    {
        #region Members
        // Static:
        private const int MaximumCopyNumber = 5;
        private const int MaxAlleleNumber = 3;
        private const int DefaultAlleleDensityThreshold = 1000;
        private const double MaxQscore = 60.0;
        private const int DefaultPerSegmentAlleleMaxCounts = 100;
        private const int DefaultAlleleCountThreshold = 4;
        private const int MaxNumOffspringGenotypes = 500;
        private const double DeNovoRate = 0.00001;
        private const int MinimumCallSize = 1000;
        // QualityFilterThreshold based on 
        public int QualityFilterThreshold { get; } = 7;
        public int DeNovoQualityFilterThreshold { get; } = 20;

        #endregion

        internal int CallVariantsInPedigree(List<string> variantFrequencyFiles, List<string> segmentFiles, List<string> outVcfFiles, string ploidyBedPath, string referenceFolder, List<string> sampleNames, string pedigreeFile)
        {
            // load files
            // initialize data structures and classes
            int fileCounter = 0;
            Dictionary<string, PedigreeMember.Kinship> kinships = ReadPedigreeFile(pedigreeFile);
            List<PedigreeMember> pedigreeMembers = new List<PedigreeMember>();
            foreach (string sampleName in sampleNames)
            {
                var pedigreeMember = SetPedigreeMember(variantFrequencyFiles, segmentFiles, ploidyBedPath, sampleName, fileCounter);
                pedigreeMember.Kin = kinships[pedigreeMember.Name] == PedigreeMember.Kinship.Parent ?
                    PedigreeMember.Kinship.Parent : PedigreeMember.Kinship.Offspring;
                if (kinships[pedigreeMember.Name] == PedigreeMember.Kinship.Proband)
                    pedigreeMember.Kin = PedigreeMember.Kinship.Proband;
                pedigreeMembers.Add(pedigreeMember);
                fileCounter++;
            }

            var numberOfSegments = pedigreeMembers.First().Segments.Count;
            List<GenomicInterval> segmentIntervals = GetParallelIntevals(numberOfSegments, Environment.ProcessorCount);

            List<PedigreeMember> parents = GetParents(pedigreeMembers);
            List<PedigreeMember> offsprings = GetChildren(pedigreeMembers);
            double[][] transitionMatrix = GetTransitionMatrix(MaximumCopyNumber);
            List<Tuple<int, int>> parentalGenotypes = GenerateParentalGenotypes(MaximumCopyNumber);
            List<List<Tuple<int, int>>> offspringsGenotypes = new List<List<Tuple<int, int>>>(Convert.ToInt32(Math.Pow(parentalGenotypes.Count, offsprings.Count)));
            GenerateOffspringGenotypes(offspringsGenotypes, parentalGenotypes, offsprings.Count, new List<Tuple<int, int>>());
            if (offspringsGenotypes.Count > MaxNumOffspringGenotypes)
            {
                offspringsGenotypes.Shuffle();
                offspringsGenotypes = offspringsGenotypes.Take(MaxNumOffspringGenotypes).ToList();
            }

            foreach (PedigreeMember parent in parents)
                parent.CnModel = new CopyNumberModel(MaximumCopyNumber, parent.MeanCoverage / 2.0, parent.MeanMafCoverage / 2.0, parent.Variance, parent.MafVariance, parent.MaxCoverage);
            foreach (PedigreeMember offspring in offsprings)
                offspring.CnModel = new CopyNumberModel(MaximumCopyNumber, offspring.MeanCoverage / 2.0, offspring.MeanMafCoverage / 2.0, offspring.Variance, offspring.MafVariance, offspring.MaxCoverage);

            Parallel.ForEach(
                segmentIntervals,
                new ParallelOptions
                {
                    MaxDegreeOfParallelism = Environment.ProcessorCount,
                    TaskScheduler = TaskScheduler.Default
                },
                interval =>
                {
                    Console.WriteLine($"{DateTime.Now} Launching SPW task for segment {interval.Start} - {interval.End}");
                    var segmentIndex = 0;
                    while (segmentIndex < numberOfSegments)
                    {
                        if (segmentIndex >= interval.Start && segmentIndex <= interval.End)
                        {
                            var alleleCounts = pedigreeMembers.Select(x => x.Segments[segmentIndex].Alleles.Counts.Count);
                            var enumerable = alleleCounts as int[] ?? alleleCounts.ToArray();
                            var alleleDensity = enumerable.Average() / pedigreeMembers.First().Segments[segmentIndex].Length;
                            var useCnLikelihood = enumerable.Select(x => x > DefaultAlleleCountThreshold).Any(c => c == false) &&
                                alleleDensity < DefaultAlleleDensityThreshold && enumerable.Average() < DefaultPerSegmentAlleleMaxCounts;
                            CopyNumberDistribution copyNumberLikelihoods = useCnLikelihood ?
                            MaximalCnLikelihoodWithPedigreeInfo(parents, offsprings, segmentIndex, transitionMatrix, offspringsGenotypes) :
                            MaximalGtLikelihoodWithPedigreeInfo(parents, offsprings, segmentIndex, parentalGenotypes, offspringsGenotypes);
                            EstimateQScoresWithPedigreeInfo(parents, offsprings, segmentIndex, copyNumberLikelihoods);
                        }

                        segmentIndex++;
                    }
                    Console.WriteLine($"{DateTime.Now} Finished SPW task for segment {interval.Start} - {interval.End}");
                });

            int pedigreeMemberIndex = 0;
            foreach (var pedigreeMember in pedigreeMembers)
            {
                CanvasSegment.MergeSegments(ref pedigreeMember.Segments, MinimumCallSize);
                CanvasSegment.WriteSegments(outVcfFiles[pedigreeMemberIndex], pedigreeMember.Segments,
                    pedigreeMember.MeanCoverage, referenceFolder, pedigreeMember.Name, null, null,
                    QualityFilterThreshold, DeNovoQualityFilterThreshold);
                pedigreeMemberIndex++;
            }
            return 0;
        }

        internal int CallVariants(List<string> variantFrequencyFiles, List<string> segmentFiles, List<string> outVcfFiles, string ploidyBedPath, string referenceFolder, List<string> sampleNames, string pedigreeFile)
        {
            // load files
            // initialize data structures and classes
            int fileCounter = 0;
            List<PedigreeMember> pedigreeMembers = new List<PedigreeMember>();
            foreach (string sampleName in sampleNames)
            {
                var pedigreeMember = SetPedigreeMember(variantFrequencyFiles, segmentFiles, ploidyBedPath, sampleName, fileCounter);
                pedigreeMembers.Add(pedigreeMember);
                fileCounter++;
            }

            var numberOfSegments = pedigreeMembers.First().Segments.Count;
            List<GenomicInterval> segmentIntervals = GetParallelIntevals(numberOfSegments, Environment.ProcessorCount);
            var genotypes = GenerateGenotypeCombinations(MaximumCopyNumber, MaxAlleleNumber);
            var copyNumberCombinations = GenerateCopyNumberCombinations(MaximumCopyNumber, MaxAlleleNumber);


            foreach (PedigreeMember pedigreeMember in pedigreeMembers)
                pedigreeMember.CnModel = new CopyNumberModel(MaximumCopyNumber, pedigreeMember.MeanCoverage / 2.0, pedigreeMember.MeanMafCoverage / 2.0,
                    pedigreeMember.Variance, pedigreeMember.MafVariance, pedigreeMember.MaxCoverage);

            Parallel.ForEach(
                segmentIntervals,
                new ParallelOptions
                {
                    MaxDegreeOfParallelism = Environment.ProcessorCount,
                    TaskScheduler = TaskScheduler.Default
                },
                interval =>
                {
                    Console.WriteLine($"{DateTime.Now} Launching SPW task for segment {interval.Start} - {interval.End}");
                    var segmentIndex = 0;
                    while (segmentIndex < numberOfSegments)
                    {
                        if (segmentIndex >= interval.Start && segmentIndex <= interval.End)
                        {
                            var alleleCounts = pedigreeMembers.Select(x => x.Segments[segmentIndex].Alleles.Counts.Count);
                            var enumerable = alleleCounts as int[] ?? alleleCounts.ToArray();
                            var alleleDensity = enumerable.Average() / pedigreeMembers.First().Segments[segmentIndex].Length;
                            var useCnLikelihood = enumerable.Select(x => x > DefaultAlleleCountThreshold).Any(c => c == false) &&
                                alleleDensity < DefaultAlleleDensityThreshold && enumerable.Average() < DefaultPerSegmentAlleleMaxCounts;
                            double[][] copyNumberLikelihoods = useCnLikelihood ?
                            MaximalCnLikelihoodNoPedigreeInfo(pedigreeMembers, segmentIndex, copyNumberCombinations) :
                            MaximalGtLikelihoodNoPedigreeInfo(pedigreeMembers, segmentIndex, genotypes);
                            EstimateQScoresNoPedigreeInfo(pedigreeMembers, segmentIndex, copyNumberLikelihoods);
                        }
                        segmentIndex++;
                    }
                    Console.WriteLine($"{DateTime.Now} Finished SPW task for segment {interval.Start} - {interval.End}");
                });

            int nSegments = pedigreeMembers.First().Segments.Count;
            var copyNumbers = new List<List<int>>(nSegments);
            var qscores = new List<double>(nSegments);
            foreach (int segmentIndex in Enumerable.Range(0, nSegments))
            {
                copyNumbers.Add(new List<int>());
                qscores.Add(0);
                foreach (var pedigreeMember in pedigreeMembers)
                {
                    copyNumbers[segmentIndex].Add(pedigreeMember.Segments[segmentIndex].CopyNumber);
                    qscores[segmentIndex] += pedigreeMember.Segments[segmentIndex].QScore;
                }
                qscores[segmentIndex] /= pedigreeMembers.Count;
            }

            int pedigreeMemberIndex = 0;
            if ((copyNumbers == null && qscores != null) || (copyNumbers != null & qscores == null))
                throw new ArgumentException("Both copyNumbers and qscores arguments must be specified.");
            if (copyNumbers != null && copyNumbers.Count != pedigreeMembers.First().Segments.Count)
                throw new ArgumentException("Length of copyNumbers list should be equal to the number of segments.");
            if (qscores != null && qscores.Count != pedigreeMembers.First().Segments.Count)
                throw new ArgumentException("Length of qscores list should be equal to the number of segments.");

            foreach (var pedigreeMember in pedigreeMembers)
            {
                CanvasSegment.MergeSegments(ref pedigreeMember.Segments, MinimumCallSize, copyNumbers: copyNumbers, qscores: qscores);
                CanvasSegment.WriteSegments(outVcfFiles[pedigreeMemberIndex], pedigreeMember.Segments,
                    pedigreeMember.MeanCoverage, referenceFolder, pedigreeMember.Name, null, null,
                    QualityFilterThreshold, DeNovoQualityFilterThreshold);
                pedigreeMemberIndex++;
            }
            return 0;
        }

        private static PedigreeMember SetPedigreeMember(List<string> variantFrequencyFiles, List<string> segmentFiles, string ploidyBedPath,
            string sampleName, int fileCounter)
        {
            PedigreeMember pedigreeMember = new PedigreeMember();
            pedigreeMember.Name = sampleName;
            pedigreeMember.Segments = CanvasSegment.ReadSegments(segmentFiles[fileCounter]);
            pedigreeMember.MeanMafCoverage = CanvasIO.LoadFrequencies(variantFrequencyFiles[fileCounter],
                pedigreeMember.Segments);
            foreach (CanvasSegment segment in pedigreeMember.Segments)
                if (segment.Alleles.Counts.Count > DefaultAlleleCountThreshold)
                    segment.Alleles.SetMedianCounts();
            pedigreeMember.Variance =
                Math.Pow(Utilities.StandardDeviation(pedigreeMember.Segments.Select(x => x.MedianCount).ToArray()), 2);
            pedigreeMember.MafVariance =
                Math.Pow(
                    Utilities.StandardDeviation(
                        pedigreeMember.Segments.Where(x => x.Alleles.TotalCoverage.Count > 0)
                            .Select(x => x.Alleles.TotalCoverage.Average())
                            .ToArray()), 2);
            pedigreeMember.MeanCoverage = pedigreeMember.Segments.Select(x => x.MedianCount).Average();
            pedigreeMember.MaxCoverage = Convert.ToInt32(pedigreeMember.Segments.Select(x => x.MedianCount).Max() + 10);
            pedigreeMember.Ploidy = PloidyInfo.LoadPloidyFromVcfFile(ploidyBedPath, pedigreeMember.Name);
            return pedigreeMember;
        }

        private void EstimateQScoresWithPedigreeInfo(List<PedigreeMember> parents, List<PedigreeMember> offsprings, int segmentIndex,
            CopyNumberDistribution copyNumberLikelihoods)
        {
            var cnStates = GetCnStates(parents, offsprings, segmentIndex);
            var names = parents.Concat(offsprings).Select(x => x.Name).ToList();
            var probands = GetProbands(offsprings);
            var singleSampleQualityScores = GetSingleSampleQualityScores(copyNumberLikelihoods, cnStates, names);

            var parent1Index = names.IndexOf(parents.First().Name);
            var parent2Index = names.IndexOf(parents.Last().Name);

            foreach (PedigreeMember proband in probands)
            {
                var probandIndex = names.IndexOf(proband.Name);
                var remainingProbandIndex = probands.Except(proband.ToSingleItemEnumerable()).Select(x => names.IndexOf(x.Name));

                if (cnStates[probandIndex] != 2 && cnStates[parent1Index] == 2 && cnStates[parent2Index] == 2 &&
                    remainingProbandIndex.All(index => cnStates[index] == 2) && singleSampleQualityScores[probandIndex] > QualityFilterThreshold)
                {
                    var deNovoQualityScore = GetConditionalDeNovoQualityScore(copyNumberLikelihoods, probandIndex,
                        cnStates[probandIndex], names[probandIndex], parent1Index, parent2Index, remainingProbandIndex.ToList());
                    if (Double.IsInfinity(deNovoQualityScore) | deNovoQualityScore > MaxQscore)
                        deNovoQualityScore = MaxQscore;
                    probands.First().Segments[segmentIndex].DQScore = deNovoQualityScore;
                }
            }

            var counter = 0;
            foreach (PedigreeMember sample in parents.Concat(offsprings))
            {
                sample.Segments[segmentIndex].QScore = singleSampleQualityScores[counter];
                if (sample.Segments[segmentIndex].QScore < QualityFilterThreshold)
                    sample.Segments[segmentIndex].Filter = $"q{QualityFilterThreshold}";
                counter++;
            }
        }

        private void EstimateQScoresNoPedigreeInfo(List<PedigreeMember> samples, int segmentIndex, double[][] copyNumberLikelihoods)
        {
            var cnStates = samples.Select(x => Math.Min(x.Segments[segmentIndex].CopyNumber, MaximumCopyNumber - 1)).ToList();
            int counter = 0;
            foreach (PedigreeMember sample in samples)
            {
                double normalizationConstant = copyNumberLikelihoods[counter].Sum();
                var qscore = -10.0 * Math.Log10((normalizationConstant - copyNumberLikelihoods[counter][cnStates[counter]]) / normalizationConstant);
                if (Double.IsInfinity(qscore) | qscore > MaxQscore)
                    qscore = MaxQscore;
                sample.Segments[segmentIndex].QScore = qscore;
                if (sample.Segments[segmentIndex].QScore < QualityFilterThreshold)
                    sample.Segments[segmentIndex].Filter = $"q{QualityFilterThreshold}";
                counter++;
            }
        }


        private static List<int> GetCnStates(IEnumerable<PedigreeMember> parents, IEnumerable<PedigreeMember> offsprings, int segmentIndex)
        {
            return parents.Concat(offsprings).Select(x => Math.Min(x.Segments[segmentIndex].CopyNumber, MaximumCopyNumber - 1)).ToList();
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
        public CopyNumberDistribution MaximalCnLikelihoodWithPedigreeInfo(List<PedigreeMember> parents, List<PedigreeMember> children, int segmentPosition, double[][] transitionMatrix, List<List<Tuple<int, int>>> offspringsGenotypes)
        {
            double maximalLikelihood;
            InitializeLikelihood(out maximalLikelihood, segmentPosition, parents, children);
            var parent1Likelihood = parents.First().CnModel.GetCnLikelihood(Math.Min(parents.First().GetCoverage(segmentPosition), parents.First().MeanCoverage * 3.0));
            var parent2Likelihood = parents.Last().CnModel.GetCnLikelihood(Math.Min(parents.Last().GetCoverage(segmentPosition), parents.Last().MeanCoverage * 3.0));

            if (parent1Likelihood.Count != parent2Likelihood.Count)
                throw new ArgumentException("Both parents should have the same number of CN states");
            int nCopies = MaximumCopyNumber;
            List<string> names = parents.Select(x => x.Name).Union(children.Select(x => x.Name)).ToList();
            var density = new CopyNumberDistribution(nCopies, names);

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
                            parents.Last().Segments[segmentPosition].CopyNumber = cn2;
                            counter = 0;
                            foreach (PedigreeMember child in children)
                            {
                                child.Segments[segmentPosition].CopyNumber = offspringGtStates[counter].Item1 + offspringGtStates[counter].Item2;
                                counter++;
                            }
                        }
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
        public CopyNumberDistribution MaximalGtLikelihoodWithPedigreeInfo(List<PedigreeMember> parents, List<PedigreeMember> children, int segmentPosition, List<Tuple<int, int>> parentalGenotypes, List<List<Tuple<int, int>>> offspringsGenotypes)
        {
            double maximalLikelihood;
            InitializeLikelihood(out maximalLikelihood, segmentPosition, parents, children);
            var parent1Likelihood = parents.First().CnModel.GetMedianGtLikelihood(parents.First().GetMedianAlleleCounts(segmentPosition));
            var parent2Likelihood = parents.Last().CnModel.GetMedianGtLikelihood(parents.Last().GetMedianAlleleCounts(segmentPosition));
            int nCopies = MaximumCopyNumber;
            List<string> names = parents.Select(x => x.Name).Union(children.Select(x => x.Name)).ToList();
            var density = new CopyNumberDistribution(nCopies, names);

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
                            currentLikelihood *= GetTransitionProbability(parent1GtStates.Item1, parent1GtStates.Item2, offspringGtStates[counter].Item1, offspringGtStates[counter].Item2) *
                                                 GetTransitionProbability(parent2GtStates.Item1, parent2GtStates.Item2, offspringGtStates[counter].Item1, offspringGtStates[counter].Item2) *
                                                 child.CnModel.GetMedianGtLikelihood(child.GetMedianAlleleCounts(segmentPosition))[offspringGtStates[counter].Item1][offspringGtStates[counter].Item2];
                            counter++;
                        }

                        int[] copyNumberIndices = { parent1GtStates.Item1 + parent1GtStates.Item2, parent2GtStates.Item1 + parent2GtStates.Item2 };
                        var index = copyNumberIndices.Concat(offspringGtStates.Select(x => x.Item1 + x.Item2)).ToArray();
                        density.SetJointProbability(Math.Max(currentLikelihood, density.GetJointProbability(index)), index);

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
                    }
                }
            }
            return density;
        }

        /// <summary>
        /// Calculates maximal likelihood for segments without SNV allele ratios. Updated CanvasSegment CopyNumber only. 
        /// </summary>
        /// <param name="parents"></param>
        /// <param name="children"></param>
        /// <param name="segmentPosition"></param>
        /// <param name="transitionMatrix"></param>
        public double[][] MaximalCnLikelihoodNoPedigreeInfo(List<PedigreeMember> samples, int segmentPosition, List<List<int>> copyNumberCombinations)
        {
            int defaultCn = 2;
            double maximalLikelihood = 0;
            foreach (PedigreeMember sample in samples)
                sample.Segments[segmentPosition].CopyNumber = defaultCn;
            int nCopies = MaximumCopyNumber;
            List<string> names = samples.Select(x => x.Name).ToList();
            var totalLikelihoods = new List<double>();
            foreach (var copyNumberCombination in copyNumberCombinations)
            {
                double totalLikelihood = 0;
                foreach (PedigreeMember sample in samples)
                {
                    maximalLikelihood = 0;
                    foreach (var copyNumber in copyNumberCombination)
                    {
                        var currentLikelihood = sample.CnModel.GetCnLikelihood(Math.Min(sample.GetCoverage(segmentPosition), sample.MeanCoverage * 3.0))[copyNumber];
                        currentLikelihood = Double.IsNaN(currentLikelihood) || Double.IsInfinity(currentLikelihood) ? 0 : currentLikelihood;
                        if (currentLikelihood > maximalLikelihood)
                            maximalLikelihood = currentLikelihood;
                    }
                    totalLikelihood += maximalLikelihood;
                }
                totalLikelihoods.Add(totalLikelihood);
            }

            var bestcopyNumberCombination = copyNumberCombinations[totalLikelihoods.IndexOf(totalLikelihoods.Max())];
            int counter = 0;
            var density = new double[samples.Count][];
            foreach (PedigreeMember sample in samples)
            {
                maximalLikelihood = 0;
                density[counter] = new double[nCopies];
                counter++;
                foreach (var copyNumber in bestcopyNumberCombination)
                {
                    var currentLikelihood = sample.CnModel.GetCnLikelihood(Math.Min(sample.GetCoverage(segmentPosition), sample.MeanCoverage * 3.0))[copyNumber];
                    currentLikelihood = Double.IsNaN(currentLikelihood) || Double.IsInfinity(currentLikelihood) ? 0 : currentLikelihood;
                    if (currentLikelihood > maximalLikelihood)
                    {
                        maximalLikelihood = currentLikelihood;
                        sample.Segments[segmentPosition].CopyNumber = copyNumber;
                        density[names.FindIndex(name => name == sample.Name)][copyNumber] = maximalLikelihood;
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
        public double[][] MaximalGtLikelihoodNoPedigreeInfo(List<PedigreeMember> samples, int segmentPosition,
            List<Dictionary<int, List<Tuple<int, int>>>> genotypesets)
        {
            int defaultCn = 2;
            double maximalLikelihood = 0;
            foreach (PedigreeMember sample in samples)
                sample.Segments[segmentPosition].CopyNumber = defaultCn;
            int nCopies = MaximumCopyNumber;
            List<string> names = samples.Select(x => x.Name).ToList();
            var totalLikelihoods = new List<double>();
            foreach (var genotypeset in genotypesets)
            {
                double totalLikelihood = 0;
                foreach (PedigreeMember sample in samples)
                {
                    maximalLikelihood = 0;
                    foreach (var genotype in genotypeset)
                    {
                        var selectedGtState = 0;
                        var currentLikelihood = sample.CnModel.GetGtLikelihood(sample.GetAlleleCounts(segmentPosition),
                            genotype.Value, ref selectedGtState, sample.MaxCoverage);
                        currentLikelihood = Double.IsNaN(currentLikelihood) || Double.IsInfinity(currentLikelihood) ? 0 : currentLikelihood;
                        if (currentLikelihood > maximalLikelihood)
                            maximalLikelihood = currentLikelihood;
                    }
                    totalLikelihood += maximalLikelihood;
                }
                totalLikelihoods.Add(totalLikelihood);
            }

            var bestGenotypeset = genotypesets[totalLikelihoods.IndexOf(totalLikelihoods.Max())];
            int counter = 0;
            var density = new double[samples.Count][];

            foreach (PedigreeMember sample in samples)
            {
                maximalLikelihood = 0;
                density[counter] = new double[nCopies];
                counter++;
                foreach (var genotype in bestGenotypeset)
                {
                    var selectedGtState = 0;
                    var currentLikelihood = sample.CnModel.GetGtLikelihood(sample.GetAlleleCounts(segmentPosition),
                        genotype.Value, ref selectedGtState, sample.MaxCoverage);
                    int copyNumber = genotype.Key;
                    currentLikelihood = Double.IsNaN(currentLikelihood) || Double.IsInfinity(currentLikelihood)? 0 : currentLikelihood;
                    density[names.FindIndex(name => name == sample.Name)][copyNumber] = Math.Max(currentLikelihood,
                        density[names.FindIndex(name => name == sample.Name)][copyNumber]);

                    if (currentLikelihood > maximalLikelihood)
                    {
                        maximalLikelihood = currentLikelihood;
                        sample.Segments[segmentPosition].CopyNumber = copyNumber;
                        sample.Segments[segmentPosition].MajorChromosomeCount = Math.Max(genotype.Value[selectedGtState].Item1,
                            genotype.Value[selectedGtState].Item2);
                    }
                }
            }
            return density;
        }


        private static void InitializeLikelihood(out double maximalLikelihood, int segmentPosition, List<PedigreeMember> parents, List<PedigreeMember> children)
        {
            maximalLikelihood = 0;
            int defaultCn = 2;

            parents.First().Segments[segmentPosition].CopyNumber = defaultCn;
            parents.Last().Segments[segmentPosition].CopyNumber = defaultCn;
            foreach (PedigreeMember child in children)
                child.Segments[segmentPosition].CopyNumber = defaultCn;
        }

        /// <summary>
        /// Generate all possible copy number combinations with the maximal number of copy numbers per segment set to maxAlleleNumber.
        /// </summary>
        /// <param name="numberOfCnStates"></param>
        /// <param name="maxAlleleNumber"></param>
        /// <returns></returns>
        public static List<List<int>> GenerateCopyNumberCombinations(int numberOfCnStates, int maxAlleleNumber)
        {
            if (numberOfCnStates <= 0)
                throw new ArgumentOutOfRangeException(nameof(numberOfCnStates));
            var cnStates = Enumerable.Range(0, numberOfCnStates).ToList();
            var allCombinations = new List<List<int>>();
            for (int currentAlleleNumber = 1; currentAlleleNumber <= maxAlleleNumber; currentAlleleNumber++)
            {
                var permutations = new Combinations<int>(cnStates, currentAlleleNumber);
                var list = permutations.Select(x => x.ToList()).ToList();
                allCombinations.AddRange(list);
            }
            return allCombinations;
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

        /// <summary>
        /// Generate all possible copy number genotype combinations with the maximal number of alleles per segment set to maxAlleleNumber.
        /// </summary>
        /// <param name="numberOfCnStates"></param>
        /// <param name="maxAlleleNumber"></param>
        /// <returns> </returns>
        public List<Dictionary<int, List<Tuple<int, int>>>> GenerateGenotypeCombinations(int numberOfCnStates, int maxAlleleNumber)
        {
            Dictionary<int, List<Tuple<int, int>>> genotypes = new Dictionary<int, List<Tuple<int, int>>>();
            for (int cn = 0; cn < numberOfCnStates; cn++)
            {
                genotypes[cn] = new List<Tuple<int, int>>();
                for (int gt = 0; gt <= cn; gt++)
                {
                    genotypes[cn].Add(new Tuple<int, int>(gt, cn - gt));
                }
            }

            var cnStates = Enumerable.Range(0, numberOfCnStates).ToList();
            var allCombinations = new List<List<int>>();
            for (int currentAlleleNumber = 1; currentAlleleNumber <= maxAlleleNumber; currentAlleleNumber++)
            {
                var permutations = new Combinations<int>(cnStates, currentAlleleNumber);
                var list = permutations.Select(x => x.ToList()).ToList();
                allCombinations.AddRange(list);
            }
            List < Dictionary < int, List < Tuple<int, int>>>> genotypeCombinations = new List<Dictionary<int, List<Tuple<int, int>>>>();
            foreach (List<int> combination in allCombinations)
            {
                var tmpDictionary = new Dictionary<int, List<Tuple<int, int>>>();
                foreach (int index in combination)
                {
                    int key = genotypes.Keys.ToList()[index];
                    tmpDictionary.Add(key, genotypes[index]);
                }
                genotypeCombinations.Add(tmpDictionary);
            }
            return genotypeCombinations;
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
            int step = nSegments / nCores;
            intevals.Add(new GenomicInterval(0, step));
            int cumSum = step + 1;
            while (cumSum + step + 1 < nSegments - 1)
            {
                intevals.Add(new GenomicInterval(cumSum, cumSum + step));
                cumSum += step + 1;
            }
            intevals.Add(new GenomicInterval(cumSum, nSegments - 1));
            return intevals;
        }

        public double GetTransitionProbability(int gt1Parent, int gt2Parent, int gt1Offspring, int gt2Offspring)
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
                var qscore = -10.0 * Math.Log10((normalizationConstant - cnMarginalProbabilities[cnStates[index]]) / normalizationConstant);
                if (Double.IsInfinity(qscore) | qscore > MaxQscore)
                    qscore = MaxQscore;
                singleSampleQualityScores.Add(qscore);
            }
            return singleSampleQualityScores;
        }

        public double GetDeNovoQualityScore(List<PedigreeMember> parents, CopyNumberDistribution density, string sampleName, int sampleValue, double sampleProbability)
        {
            int nSamples = density.Count;
            int diploidState = 2;
            var probandMarginalProbabilities = density.GetMarginalProbability(nSamples, MaximumCopyNumber, sampleName);
            var normalization = probandMarginalProbabilities[sampleValue] + probandMarginalProbabilities[diploidState];
            var probandMarginalAlt = probandMarginalProbabilities[sampleValue] / normalization;
            //density.SetConditionalProbability(density.Count, MaximumCopyNumber, sampleName, sampleValue, probandMarginalProbabilities[sampleValue]);

            var parentNames = parents.Select(x => x.Name).ToList();
            var firstParentMarginalProbabilities = density.GetMarginalProbability(nSamples, MaximumCopyNumber, parentNames.First());
            var secondParentMarginalProbabilities = density.GetMarginalProbability(nSamples, MaximumCopyNumber, parentNames.Last());
            normalization = firstParentMarginalProbabilities[sampleValue] + firstParentMarginalProbabilities[diploidState];
            var firstParentMarginalAlt = Math.Min(Math.Max(firstParentMarginalProbabilities[sampleValue] / normalization, 0.001), 0.999);
            normalization = secondParentMarginalProbabilities[sampleValue] + secondParentMarginalProbabilities[diploidState];
            var secondParentMarginalAlt = Math.Min(Math.Max(secondParentMarginalProbabilities[sampleValue] / normalization, 0.001), 0.999);

            normalization = (1 - firstParentMarginalAlt) * secondParentMarginalAlt + firstParentMarginalAlt * secondParentMarginalAlt + (1 - firstParentMarginalAlt) * (1 - firstParentMarginalAlt) +
            (1 - secondParentMarginalAlt) * firstParentMarginalAlt;
            var diploidProbability = (1 - firstParentMarginalAlt) * (1 - secondParentMarginalAlt) / normalization;
            var denovoProbability = diploidProbability * probandMarginalAlt;
            var qscore = -10.0 * Math.Log10(1 - denovoProbability);
            return qscore;
        }

        public double GetConditionalDeNovoQualityScore(CopyNumberDistribution density, int probandIndex,
                    int probandCopyNumber, string probandName, int parent1Index, int parent2Index, List<int> remainingProbandIndex)
        {

            double numerator = 0.0;
            double denominator = 0.0;
            const int diploidState = 2;
            int nSamples = density.Count;
            var probandMarginalProbabilities = density.GetMarginalProbability(nSamples, MaximumCopyNumber, probandName);
            var normalization = probandMarginalProbabilities[probandCopyNumber] + probandMarginalProbabilities[diploidState];
            var probandMarginalAlt = probandMarginalProbabilities[probandCopyNumber] / normalization;

            foreach (var copyNumberIndex in density.Indices.Where(x => x[probandIndex] == probandCopyNumber).ToArray())
            {
                if (density.GetJointProbability(copyNumberIndex.ToArray()) > 0.0)
                {
                    var holder = density.GetJointProbability(copyNumberIndex.ToArray());
                    denominator += holder;

                    if (copyNumberIndex[parent1Index] == diploidState && copyNumberIndex[parent2Index] == diploidState && remainingProbandIndex.All(index => copyNumberIndex[index] == 2))
                        numerator += holder;
                }
            }

            const double q60 = 0.000001;
            var denovoProbability = (1 - numerator / denominator) * (1 - probandMarginalAlt);
            var qscore = -10.0 * Math.Log10(Math.Max(denovoProbability, q60));
            return qscore;
        }
    }

}
