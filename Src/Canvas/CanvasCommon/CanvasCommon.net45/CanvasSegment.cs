using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Isas.SequencingFiles;

namespace CanvasCommon
{
    public class Alleles
    {
        public List<float> Frequencies = new List<float>();
        public List<int> TotalCoverage = new List<int>();
        public List<Tuple<int, int>> Counts = new List<Tuple<int, int>>();
        public Tuple<int, int> MedianCounts = new Tuple<int, int>(0, 0);
        public void SetMedianCounts()
        {
            var item1 = Utilities.Median(Counts.Select(x => Math.Max(x.Item1, x.Item2)).ToList());
            var item2 = Utilities.Median(Counts.Select(x => Math.Min(x.Item1, x.Item2)).ToList());
            MedianCounts = new Tuple<int, int>(item1, item2);
        }
    }
    /// <summary>
    /// Contains information about a genomic interval. Has functions for computing copy numbers and their likelihoods.
    /// </summary>
    public class CanvasSegment
    {
        #region Members
        public List<float> Counts;
        public int CopyNumber { get; set; }
        public int SecondBestCopyNumber { get; set; }
        public int? MajorChromosomeCount;
        public double QScore;
        public double? DQScore;
        public double? MajorChromosomeCountScore;
        public double ModelDistance;
        public double RunnerUpModelDistance;
        public bool CopyNumberSwapped;
        public bool IsHeterogeneous;
        private static readonly int NumberVariantFrequencyBins = 100;
        public string Filter = "PASS";
        public Tuple<int, int> StartConfidenceInterval; // if not null, this is a confidence interval around Start, reported in the CIPOS tag
        public Tuple<int, int> EndConfidenceInterval; // if not null, this is a confidence interval around End, reported in the CIEND tag
        public Alleles Alleles;
        public string Chr { get; private set; }
        /// <summary>
        /// bed format start position
        /// zero-based inclusive start position
        /// </summary>
        public int Begin { get; private set; }
        /// <summary>
        /// bed format end position
        /// zero-based exclusive end position (i.e. the same as one-based inclusive end position)
        /// </summary>
        public int End { get; private set; }
        #endregion
        public int Length => End - Begin;
        /// <summary>
        /// Mean of the segment's counts.
        /// </summary>
        public double MeanCount
        {
            get
            {
                double sum = this.Counts.Sum();
                return sum / this.BinCount;
            }
        }

        public double? MeanMAF
        {
            get
            {
                if (this.Alleles.Frequencies.Count <= 5)
                    return null;

                return this.Alleles.Frequencies.Select(VF => VF > 0.5 ? 1 - VF : VF).Average();
            }
        }

        /// <summary>
        /// Median of the segment's counts.
        /// </summary>
        public double MedianCount
        {
            get
            {
                var sorted = new SortedList<double>(this.Counts.Select(x => Convert.ToDouble(x)));
                return sorted.Median();
            }
        }


        /// <summary>
        /// removes flanking bins before median estimation
        /// </summary>
        public double TruncatedMedianCount(int bins2Remove)
        {
            var tmpMedian = new SortedList<double>();
            int start = Convert.ToInt32(Counts.Count * 0.25);
            int end = Counts.Count - bins2Remove;
            if (end - start > 5)
            {
                for (int index = bins2Remove; index < end; index++)
                {
                    tmpMedian.Add(Counts[index]);
                }
                return tmpMedian.Median();

            }
            var sorted = new SortedList<double>(this.Counts.Select(x => Convert.ToDouble(x)));
            return sorted.Median();
        }
        public CnvType GetCnvType(int referenceCopyNumber)
        {
            if (CopyNumber < referenceCopyNumber)
                return CnvType.Loss;
            if (CopyNumber > referenceCopyNumber)
                return CnvType.Gain;
            if (referenceCopyNumber == 2 &&
                MajorChromosomeCount.HasValue && MajorChromosomeCount == CopyNumber)
                return CnvType.LossOfHeterozygosity;
            return CnvType.Reference;
        }

        /// <summary>
        /// Merge another neighboring segment into this one.
        /// </summary>
        /// <param name="s">Segment to merge in.</param>
        public void MergeIn(CanvasSegment s)
        {
            if (s.Begin < this.Begin)
            {
                this.StartConfidenceInterval = s.StartConfidenceInterval;
                this.Begin = s.Begin;
            }
            if (s.End > this.End)
            {
                this.EndConfidenceInterval = s.EndConfidenceInterval;
                this.End = s.End;
            }
            this.Counts.AddRange(s.Counts);
            Alleles.Frequencies.AddRange(s.Alleles.Frequencies);
            Alleles.TotalCoverage.AddRange(s.Alleles.TotalCoverage);
        }

        public CanvasSegment(string chr, int begin, int end, List<float> counts)
        {
            this.Chr = chr;
            this.Begin = begin;
            this.End = end;
            this.Counts = new List<float>(counts);
            this.CopyNumber = -1;
            this.SecondBestCopyNumber = -1;
            this.Alleles = new Alleles();
        }

        public int BinCount => Counts.Count;


        /// <summary>
        /// Compute the median count from a list of segments.
        /// </summary>
        /// <param name="segments">List of segments.</param>
        /// <returns>Median of counts.</returns>
        public static double ExpectedCount(List<CanvasSegment> segments)
        {
            List<double> counts = new List<double>();

            // Get all of the counts contained within all autosomal segments.
            foreach (CanvasSegment segment in segments)
            {
                if (GenomeMetadata.SequenceMetadata.IsAutosome(segment.Chr))
                {
                    foreach (float count in segment.Counts)
                        counts.Add((double)count);
                }
            }
            double mu = Utilities.Median(counts);
            return mu;
        }

        /// <summary>
        /// Loads in data produced by CanvasPartition.exe.
        /// </summary>
        /// <param name="infile">Input file.</param>
        /// <returns>A list of segments.</returns>
        public static List<CanvasSegment> ReadSegments(string infile)
        {
            Console.WriteLine("{0} Read segments from {1}", DateTime.Now, infile);
            List<CanvasSegment> segments = new List<CanvasSegment>();

            string chr = null;
            int begin = -1;

            int previousSegmentIndex = -1;
            int previousBinStart = 0;
            int previousBinEnd = 0;
            List<float> counts = new List<float>();
            Tuple<int, int> segmentStartCI = null;
            using (GzipReader reader = new GzipReader(infile))
            {
                string row = null;

                while ((row = reader.ReadLine()) != null)
                {
                    string[] fields = row.Split('\t');

                    int currentSegmentIndex = Convert.ToInt32(fields[4]);
                    int newBinStart = Convert.ToInt32(fields[1]);
                    int newBinEnd = Convert.ToInt32(fields[2]);

                    // We've moved to a new segment
                    if (currentSegmentIndex != previousSegmentIndex)
                    {
                        // Make a segment
                        if (previousSegmentIndex != -1)
                        {
                            CanvasSegment segment = new CanvasSegment(chr, begin, previousBinEnd, counts);
                            // Prepare the confidence interval for the end of the segment that just ended, based on the size of its last bin
                            // (and, if the segments abut, based on the size of the next segment's first bin):
                            int CIEnd1 = -(previousBinEnd - previousBinStart) / 2;
                            int CIEnd2 = -CIEnd1;
                            if (previousBinEnd == newBinStart)
                            {
                                CIEnd2 = (newBinEnd - newBinStart) / 2;
                            }
                            segment.EndConfidenceInterval = new Tuple<int, int>(CIEnd1, CIEnd2);
                            segment.StartConfidenceInterval = segmentStartCI;
                            segments.Add(segment);
                            counts.Clear();

                            // Prepare the confidence interval for the start of the segment that just started, based on the size of its first
                            // bin (and, if the segments abut, based on the size of the previous segment's last bin):
                            int CIStart2 = (newBinEnd - newBinStart) / 2;
                            int CIStart1 = -CIStart2;
                            if (previousBinEnd == newBinStart)
                            {
                                CIStart1 = -(previousBinEnd - previousBinStart) / 2;
                            }
                            segmentStartCI = new Tuple<int, int>(CIStart1, CIStart2);
                        }
                        else
                        {
                            int interval = (newBinEnd - newBinStart) / 2;
                            segmentStartCI = new Tuple<int, int>(-interval, interval);
                        }
                        chr = fields[0];
                        begin = Convert.ToInt32(fields[1]);
                        previousSegmentIndex = currentSegmentIndex;
                    }
                    previousBinStart = newBinStart;
                    previousBinEnd = newBinEnd;

                    counts.Add(float.Parse(fields[3]));
                }

                if (previousSegmentIndex != -1)
                {
                    // Add the last segment
                    CanvasSegment segment = new CanvasSegment(chr, begin, previousBinEnd, counts);
                    segments.Add(segment);
                    segment.StartConfidenceInterval = segmentStartCI;
                }
            }
            Console.WriteLine("{0} Loaded {1} segments", DateTime.Now, segments.Count);
            return segments;
        }



        /// <summary>
        /// Apply quality scores.
        /// </summary>
        public static void AssignQualityScores(List<CanvasSegment> segments, QScoreMethod qscoreMethod, QualityScoreParameters qscoreParameters)
        {
            foreach (CanvasSegment segment in segments)
            {
                segment.QScore = segment.ComputeQScore(qscoreMethod, qscoreParameters);
            }
        }


        static public Dictionary<string, List<CanvasSegment>> GetSegmentsByChromosome(List<CanvasSegment> segments)
        {
            Dictionary<string, List<CanvasSegment>> segmentsByChromosome = new Dictionary<string, List<CanvasSegment>>();
            foreach (CanvasSegment segment in segments)
            {
                if (!segmentsByChromosome.ContainsKey(segment.Chr))
                {
                    segmentsByChromosome[segment.Chr] = new List<CanvasSegment>();
                }
                segmentsByChromosome[segment.Chr].Add(segment);
            }
            return segmentsByChromosome;
        }

        /// <summary>
        /// Work out how many counts we expect to see for a "typical" plot point, and exclude plotting
        /// of points which have very little data available (e.g. parts of chrY), given our coverage and bin size.
        /// Old hard-coded cutoff was 30.
        /// </summary>
        /// <returns></returns>
        private static int GetMinimumBinsForCoveragePlotPoint(List<CanvasSegment> segments, int pointLength)
        {
            long totalBins = 0;
            long totalLength = 0;
            foreach (var segment in segments)
            {
                totalBins += segment.Counts.Count;
                totalLength += segment.End - segment.Begin;
            }

            // Plot points that have at least 25% as much coverage info as we expect to see on average
            return Math.Max(1, (int)(0.25f * totalBins / (totalLength / pointLength)));
        }

        /// <summary>
        /// Generate a tabular file with information about coverage and allele frequency for each chunk of the genome.
        /// This file can be used to generate a pretty plot of coverage versus MAF.  
        /// </summary>
        public static void WriteCoveragePlotData(List<CanvasSegment> segments, double? normalDiploidCoverage, PloidyInfo referencePloidy,
            string filePath, string referenceFolder)
        {
            int pointLength = 100000;
            int minimumBinsToPlot = GetMinimumBinsForCoveragePlotPoint(segments, pointLength);

            Dictionary<string, List<CanvasSegment>> segmentsByChromosome = GetSegmentsByChromosome(segments);
            GenomeMetadata genome = new GenomeMetadata();
            genome.Deserialize(new FileLocation(Path.Combine(referenceFolder, "GenomeSize.xml")));


            List<float> counts = new List<float>();
            List<float> MAF = new List<float>();
            List<float> VF = new List<float>();
            using (FileStream stream = new FileStream(filePath, FileMode.Create, FileAccess.Write))
            using (StreamWriter writer = new StreamWriter(stream))
            {
                writer.NewLine = "\n";
                writer.Write("#Chromosome\tStart\tEnd\tCopyNumber\tMajorChromosomeCount\tMedianHits\tNormalizedCoverage\tMedianMinorAlleleFrequency\tReferencePloidy\t");
                for (int i = 0; i < NumberVariantFrequencyBins; i++) { writer.Write("VariantFrequencyBin{0}\t", i); }
                writer.WriteLine();
                foreach (GenomeMetadata.SequenceMetadata chromosome in genome.Sequences)
                {
                    if (!segmentsByChromosome.ContainsKey(chromosome.Name))
                    {
                        continue;
                    }
                    int pointStartPos = 0; // 0-based start
                    while (pointStartPos < chromosome.Length)
                    {
                        int pointEndPos = (int)Math.Min(chromosome.Length, pointStartPos + pointLength); // 1-based end
                        counts.Clear();
                        MAF.Clear();
                        VF.Clear();
                        Dictionary<string, long> CopyNumberAndChromCount = new Dictionary<string, long>();
                        Dictionary<int, long> basesByCopyNumber = new Dictionary<int, long>();
                        // Accumulate counts and MAF from the segments:
                        List<CanvasSegment> chrSegments = new List<CanvasSegment>();
                        if (segmentsByChromosome.ContainsKey(chromosome.Name)) chrSegments = segmentsByChromosome[chromosome.Name];
                        List<CanvasSegment> overlapSegments = new List<CanvasSegment>();
                        foreach (CanvasSegment segment in chrSegments)
                        {
                            if (segment.Begin > pointEndPos) continue;
                            if (segment.End < pointStartPos) continue;

                            int weight = Math.Min(segment.End, pointEndPos) - Math.Max(segment.Begin, pointStartPos);
                            string key = string.Format("{0} {1}", segment.CopyNumber, segment.MajorChromosomeCount);
                            if (!CopyNumberAndChromCount.ContainsKey(key)) CopyNumberAndChromCount[key] = 0;
                            CopyNumberAndChromCount[key] += weight;
                            if (!basesByCopyNumber.ContainsKey(segment.CopyNumber)) basesByCopyNumber[segment.CopyNumber] = 0;
                            basesByCopyNumber[segment.CopyNumber] += weight;
                            overlapSegments.Add(segment);
                        }

                        // Note the most common copy number:
                        long bestCount = 0;
                        int majorCopyNumber = 0;
                        foreach (int key in basesByCopyNumber.Keys)
                        {
                            if (basesByCopyNumber[key] > bestCount)
                            {
                                bestCount = basesByCopyNumber[key];
                                majorCopyNumber = key;
                            }
                        }

                        // Find the most common major chromosome count, for the most common copy number:
                        int? majorChromosomeCount = null;
                        bestCount = 0;
                        foreach (string key in CopyNumberAndChromCount.Keys)
                        {
                            string[] bits = key.Split();
                            if (bits[1].Length == 0) continue;
                            if (int.Parse(bits[0]) != majorCopyNumber) continue;
                            long count = CopyNumberAndChromCount[key];
                            if (count < bestCount) continue;
                            bestCount = count;
                            majorChromosomeCount = int.Parse(bits[1]);
                        }

                        // Note allele frequency and coverage info, for all overlap segments that match (more or less)
                        // the most common copy number:
                        foreach (CanvasSegment segment in overlapSegments)
                        {
                            if ((majorCopyNumber == 2 && segment.CopyNumber != 2) ||
                                (majorCopyNumber < 2 && segment.CopyNumber >= 2) ||
                                (majorCopyNumber > 2 && segment.CopyNumber <= 2))
                                continue;
                            float segLength = segment.End - segment.Begin;

                            // Add counts to the overall list:
                            int firstIndex = 0;
                            if (pointStartPos > segment.Begin)
                            {
                                firstIndex = (int)((float)segment.Counts.Count * (pointStartPos - segment.Begin) / segLength);
                            }
                            int lastIndex = segment.Counts.Count;
                            if (pointEndPos < segment.End)
                            {
                                lastIndex = (int)((float)segment.Counts.Count * (pointEndPos - segment.Begin) / segLength);
                            }
                            for (int index = firstIndex; index < lastIndex; index++) counts.Add(segment.Counts[index]);

                            // Add MAF to the overall list:
                            firstIndex = 0;
                            if (pointStartPos > segment.Begin)
                            {
                                firstIndex = (int)((float)segment.Alleles.Frequencies.Count * (pointStartPos - segment.Begin) / segLength);
                            }
                            lastIndex = segment.Alleles.Frequencies.Count;
                            if (pointEndPos < segment.End)
                            {
                                lastIndex = (int)((float)segment.Alleles.Frequencies.Count * (pointEndPos - segment.Begin) / segLength);
                            }
                            for (int index = firstIndex; index < lastIndex; index++)
                            {
                                float tempMAF = segment.Alleles.Frequencies[index];
                                VF.Add(tempMAF);
                                if (tempMAF > 0.5) tempMAF = 1 - tempMAF;
                                MAF.Add(tempMAF);
                            }
                        }

                        // Write output for this point:
                        writer.Write("{0}\t{1}\t{2}\t", chromosome.Name, pointStartPos, pointEndPos);

                        // Write counts if we have reasonable amounts of data; write MAF if we have reasonable amounts of data.
                        // (Note: Observed that for germline data on chrY we often had well under 100 counts given the new, smaller bin size)
                        if (counts.Count >= minimumBinsToPlot)
                        {
                            writer.Write("{0}\t", majorCopyNumber);
                            writer.Write("{0}\t", majorChromosomeCount);
                            counts.Sort();
                            double medianHits = counts[counts.Count / 2];
                            writer.Write("{0:F2}\t", medianHits);
                            var normalizedCount = 2 * medianHits / normalDiploidCoverage;
                            var normalizedCountString = normalizedCount.HasValue ? $"{normalizedCount:F2}" : ".";
                            writer.Write($"{normalizedCountString}\t");
                            if (MAF.Count >= 10)
                            {
                                MAF.Sort();
                                writer.Write("{0}\t", MAF[MAF.Count / 2]);
                            }
                            else
                            {
                                writer.Write("\t");
                            }
                            int refPloidy = 2;
                            if (referencePloidy != null && referencePloidy.PloidyByChromosome.ContainsKey(chromosome.Name))
                            {
                                foreach (var interval in referencePloidy.PloidyByChromosome[chromosome.Name])
                                {
                                    if (interval.Start <= pointEndPos && interval.End >= pointStartPos)
                                    {
                                        refPloidy = interval.Ploidy;
                                    }
                                }
                            }
                            writer.Write("{0}\t", refPloidy);
                            if (VF.Count >= 10)
                            {
                                // bin VF
                                float[] vfDistribution = new float[NumberVariantFrequencyBins];
                                foreach (float vf in VF)
                                {
                                    int binNumber = Math.Min(vfDistribution.Length - 1, (int)Math.Floor(vf / 0.01));
                                    vfDistribution[binNumber]++;
                                }
                                for (int i = 0; i < vfDistribution.Length; i++)
                                {
                                    vfDistribution[i] = vfDistribution[i] / (float)VF.Count * 100.0f;
                                    writer.Write("{0:F2}\t", vfDistribution[i]);
                                }
                            }
                            else
                            {
                                for (int i = 0; i < NumberVariantFrequencyBins; i++) writer.Write("\t");
                            }
                        }
                        writer.WriteLine();
                        pointStartPos += pointLength;
                    }
                }
            }
        }

        /// <summary>
        /// Return true if we are not allowed to merge two segments separated by the interval (start, end).
        /// </summary>
        static private bool IsForbiddenInterval(string chr, int start, int end,
            Dictionary<string, List<SampleGenomicBin>> excludedIntervals)
        {
            if (excludedIntervals == null) return false;
            if (!excludedIntervals.ContainsKey(chr)) return false;
            foreach (SampleGenomicBin bin in excludedIntervals[chr])
            {
                if (bin.Start >= start && bin.Start <= end) return true;
                if (bin.Stop >= start && bin.Stop <= end) return true;
                if (bin.Start > end) break;
            }
            return false;
        }

        /// <summary>
        /// Iterates through a list of segments and merges those which have the same copy number call.
        /// Also, for segments smaller than MinimumCallSize, assimilate them into the neighbor with the best 
        /// quality score.  Two consecutive segments are considered neighbors if they're on the same chromosome
        /// and the space between them doesn't overlap with any excluded intervals.
        /// </summary>
        static public void MergeSegmentsUsingExcludedIntervals(ref List<CanvasSegment> segments, int MinimumCallSize,
            Dictionary<string, List<SampleGenomicBin>> excludedIntervals)
        {
            if (!segments.Any()) return;

            // Assimilate short segments into the *best* available neighbor:
            List<CanvasSegment> mergedSegments = new List<CanvasSegment>();
            int segmentIndex = 0;
            while (segmentIndex < segments.Count)
            {
                if (segments[segmentIndex].End - segments[segmentIndex].Begin >= MinimumCallSize)
                {
                    mergedSegments.Add(segments[segmentIndex]);
                    segmentIndex++;
                    continue;
                }
                int prevIndex = -1;
                double prevQ = 0;
                // Look back for a segment:
                for (int checkIndex = segmentIndex - 1; checkIndex > 0; checkIndex--)
                {
                    // Stop, if you jump to another chromosome, or cross a forbidden interval:
                    if (segments[checkIndex].Chr != segments[segmentIndex].Chr) break;
                    if (segments[checkIndex].End - segments[checkIndex].Begin < MinimumCallSize) continue;
                    if (IsForbiddenInterval(segments[checkIndex].Chr, segments[checkIndex].End, segments[segmentIndex].Begin, excludedIntervals)) break;
                    prevIndex = checkIndex;
                    prevQ = segments[checkIndex].QScore;
                    break;
                }
                // Look forward for a segment:
                int nextIndex = -1;
                double nextQ = 0;
                for (int checkIndex = segmentIndex + 1; checkIndex < segments.Count; checkIndex++)
                {
                    if (segments[checkIndex].Chr != segments[segmentIndex].Chr) break;
                    if (segments[checkIndex].End - segments[checkIndex].Begin < MinimumCallSize) continue;
                    if (IsForbiddenInterval(segments[checkIndex].Chr, segments[segmentIndex].End, segments[checkIndex].Begin, excludedIntervals)) break;
                    nextIndex = checkIndex;
                    nextQ = segments[checkIndex].QScore;
                    break;
                }

                if (prevQ > 0 && prevQ >= nextQ)
                {
                    // segments[prevIndex] assimilates segments[prevIndex+1...segmentIndex].
                    // Assimilation of previous segments was already done, so we just need to assimilate this one:
                    segments[prevIndex].MergeIn(segments[segmentIndex]);
                    segmentIndex++;
                    continue;
                }

                if (nextQ > 0)
                {
                    // segments[nextIndex] assimilates segments[segmentIndex...nextIndex - 1]
                    for (int tempIndex = segmentIndex; tempIndex < nextIndex; tempIndex++)
                    {
                        segments[nextIndex].MergeIn(segments[tempIndex]);
                    }
                    segmentIndex = nextIndex;
                    continue;
                }

                mergedSegments.Add(segments[segmentIndex]);
                segmentIndex++;
            }
            segments = mergedSegments;

            // Now, merge together adjacent segments with same calls!
            mergedSegments = new List<CanvasSegment>();
            CanvasSegment lastSegment = segments[0];
            mergedSegments.Add(lastSegment);
            segmentIndex = 1;
            while (segmentIndex < segments.Count)
            {
                // Assimilate an adjacent segment with the same copy number call and heterogeneity flag:
                if (lastSegment.CopyNumber == segments[segmentIndex].CopyNumber && lastSegment.Chr == segments[segmentIndex].Chr &&
                    !IsForbiddenInterval(lastSegment.Chr, lastSegment.End, segments[segmentIndex].Begin, excludedIntervals) &&
                    lastSegment.IsHeterogeneous == segments[segmentIndex].IsHeterogeneous)
                {
                    lastSegment.MergeIn(segments[segmentIndex]);
                    segmentIndex++;
                    continue;
                }
                lastSegment = segments[segmentIndex];
                mergedSegments.Add(segments[segmentIndex]);
                segmentIndex++;
            }
            segments = mergedSegments;
        }

        /// <summary>
        /// Iterates through a list of segments and merges those which have the same copy number call.
        /// For multisample workflow a 2D list of regions x samples is provided to test for identity of CN
        /// calls across all samples.
        /// Also, for segments smaller than MinimumCallSize, assimilate them into the neighbor with the best 
        /// quality score.  Two consecutive segments are considered neighbors if they're on the same chromosome
        /// and the space between them is not too large.
        /// </summary>
        public static void MergeSegments(ref List<CanvasSegment> segments, int minimumCallSize = 0, int maximumMergeSpan = 10000,
            List<List<int>> copyNumbers = null, List<double> qscores = null)
        {
            if (!segments.Any()) return;
            var newCopyNumbers = new List<List<int>>();

            // Assimilate short segments into the *best* available neighbor:
            List<CanvasSegment> mergedSegments = new List<CanvasSegment>();
            int segmentIndex = 0;
            while (segmentIndex < segments.Count)
            {
                if (segments[segmentIndex].End - segments[segmentIndex].Begin >= minimumCallSize)
                {
                    mergedSegments.Add(segments[segmentIndex]);
                    if (copyNumbers != null)
                        newCopyNumbers.Add(copyNumbers[segmentIndex]);
                    segmentIndex++;
                    continue;
                }
                int prevIndex = -1;
                double prevQ = -1;
                // Look back for a segment:
                for (int checkIndex = segmentIndex - 1; checkIndex >= 0; checkIndex--)
                {
                    // Stop, if you jump to another chromosome, or cross a forbidden interval:
                    if (segments[checkIndex].Chr != segments[segmentIndex].Chr) break;
                    if (segments[checkIndex].End - segments[checkIndex].Begin < minimumCallSize) continue;
                    if (segments[segmentIndex].Begin - segments[checkIndex].End > maximumMergeSpan) break;
                    prevIndex = checkIndex;
                    prevQ = qscores?[checkIndex] ?? segments[checkIndex].QScore;
                    break;
                }
                // Look forward for a segment:
                int nextIndex = -1;
                double nextQ = -1;
                for (int checkIndex = segmentIndex + 1; checkIndex < segments.Count; checkIndex++)
                {
                    if (segments[checkIndex].Chr != segments[segmentIndex].Chr) break;
                    if (segments[checkIndex].End - segments[checkIndex].Begin < minimumCallSize) continue;
                    if (segments[checkIndex].Begin - segments[segmentIndex].End > maximumMergeSpan) continue;
                    nextIndex = checkIndex;
                    nextQ = qscores?[checkIndex] ?? segments[checkIndex].QScore;
                    break;
                }

                if (prevQ >= 0 && prevQ >= nextQ)
                {
                    // segments[prevIndex] assimilates segments[prevIndex+1...segmentIndex].
                    // Assimilation of previous segments was already done, so we just need to assimilate this one:
                    segments[prevIndex].MergeIn(segments[segmentIndex]);
                    segmentIndex++;
                    continue;
                }

                if (nextQ >= 0)
                {
                    // segments[nextIndex] assimilates segments[segmentIndex...nextIndex - 1]
                    for (int tempIndex = segmentIndex; tempIndex < nextIndex; tempIndex++)
                    {
                        segments[nextIndex].MergeIn(segments[tempIndex]);
                    }
                    segmentIndex = nextIndex;
                    continue;
                }
                if (copyNumbers != null)
                    newCopyNumbers.Add(copyNumbers[segmentIndex]);
                mergedSegments.Add(segments[segmentIndex]);
                segmentIndex++;
            }
            segments = mergedSegments;
            if (copyNumbers != null && newCopyNumbers.Count != segments.Count)
                throw new ArgumentException("Length of copyNumbers list should equal the number of segments.");

            // Now, merge together adjacent segments with same calls!
            mergedSegments = new List<CanvasSegment>();
            var lastSegment = segments[0];
            var lastSegmentIndex = 0;
            mergedSegments.Add(lastSegment);
            segmentIndex = 1;
            while (segmentIndex < segments.Count)
            {
                // Assimilate an adjacent segment with the same copy number call and heterogeneity flag:
                bool mergeSegments = copyNumbers == null ? lastSegment.CopyNumber == segments[segmentIndex].CopyNumber &&
                                    lastSegment.Chr == segments[segmentIndex].Chr &&
                                    segments[segmentIndex].Begin - lastSegment.End < maximumMergeSpan &&
                                    lastSegment.IsHeterogeneous == segments[segmentIndex].IsHeterogeneous :
                                    newCopyNumbers[lastSegmentIndex].SequenceEqual(newCopyNumbers[segmentIndex]) &&
                                    lastSegment.Chr == segments[segmentIndex].Chr &&
                                    segments[segmentIndex].Begin - lastSegment.End < maximumMergeSpan;

                if (mergeSegments)
                {
                    lastSegment.MergeIn(segments[segmentIndex]);
                    segmentIndex++;
                    continue;
                }
                lastSegment = segments[segmentIndex];
                lastSegmentIndex = segmentIndex;
                mergedSegments.Add(segments[segmentIndex]);
                segmentIndex++;
            }
            segments = mergedSegments;
        }

        /// <summary>
        /// Computes QScore using one of the available methods
        /// </summary>
        public enum QScoreMethod { BinCountLinearFit, GeneralizedLinearFit, Logistic, LogisticGermline };
        public int ComputeQScore(QScoreMethod qscoreMethod, QualityScoreParameters qscoreParameters)
        {
            double score;
            int qscore;
            switch (qscoreMethod)
            {
                case QScoreMethod.LogisticGermline:
                    // Logistic model using a new selection of features.  Gives ROC curve area 0.921
                    score = qscoreParameters.LogisticGermlineIntercept;
                    score += GetQScorePredictor(QScorePredictor.LogBinCount) * qscoreParameters.LogisticGermlineLogBinCount;
                    score += GetQScorePredictor(QScorePredictor.ModelDistance) * qscoreParameters.LogisticGermlineModelDistance;
                    score += GetQScorePredictor(QScorePredictor.DistanceRatio) * qscoreParameters.LogisticGermlineDistanceRatio;
                    score = Math.Exp(score);
                    score = score / (score + 1);
                    // Transform probability into a q-score:
                    qscore = (int)(Math.Round(-10 * Math.Log10(1 - score)));
                    qscore = Math.Min(40, qscore);
                    qscore = Math.Max(2, qscore);
                    return qscore;
                case QScoreMethod.Logistic:
                    // Logistic model using a new selection of features.  Gives ROC curve area 0.8289
                    score = qscoreParameters.LogisticIntercept;
                    score += GetQScorePredictor(QScorePredictor.LogBinCount) * qscoreParameters.LogisticLogBinCount;
                    score += GetQScorePredictor(QScorePredictor.ModelDistance) * qscoreParameters.LogisticModelDistance;
                    score += GetQScorePredictor(QScorePredictor.DistanceRatio) * qscoreParameters.LogisticDistanceRatio;
                    score += GetQScorePredictor(QScorePredictor.BinCountAmpDistance);
                    score = Math.Exp(score);
                    score = score / (score + 1);
                    // Transform probability into a q-score:
                    qscore = (int)Math.Round(-10 * Math.Log10(1 - score));
                    qscore = Math.Min(60, qscore);
                    qscore = Math.Max(2, qscore);
                    return qscore;
                case QScoreMethod.BinCountLinearFit:
                    if (this.BinCount >= 100)
                        return 61;
                    else
                        return (int)Math.Round(-10 * Math.Log10(1 - 1 / (1 + Math.Exp(0.5532 - this.BinCount * 0.147))), 0, MidpointRounding.AwayFromZero);
                case QScoreMethod.GeneralizedLinearFit: // Generalized linear fit with linear transformation to QScore
                    double linearFit = qscoreParameters.GeneralizedLinearFitIntercept;
                    linearFit += qscoreParameters.GeneralizedLinearFitLogBinCount *
                                 GetQScorePredictor(QScorePredictor.LogBinCount);
                    linearFit += qscoreParameters.GeneralizedLinearFitModelDistance *
                                 GetQScorePredictor(QScorePredictor.ModelDistance);
                    linearFit += qscoreParameters.GeneralizedLinearFitMajorChromosomeCount *
                                 GetQScorePredictor(QScorePredictor.MajorChromosomeCount);
                    linearFit += qscoreParameters.GeneralizedLinearFitMafMean *
                                 GetQScorePredictor(QScorePredictor.MafMean);
                    linearFit += qscoreParameters.GeneralizedLinearFitLogMafCv * GetQScorePredictor(QScorePredictor.LogMafCv);
                    linearFit += GetQScorePredictor(QScorePredictor.BinCountAmpDistance);
                    score = -11.9 - 11.4 * linearFit; // Scaling to achieve 2 <= qscore <= 61
                    score = Math.Max(2, score);
                    score = Math.Min(61, score);
                    return (int)Math.Round(score, 0, MidpointRounding.AwayFromZero);
                default:
                    throw new Exception("Unhandled qscore method");
            }
        }

        /// <summary>
        /// Computes QScore predictor
        /// </summary>
        public enum QScorePredictor
        {
            BinCount, LogBinCount, BinCountAmpDistance, BinMean, BinCv, MafCount, MafMean, MafCv, LogMafCv, ModelDistance,
            RunnerUpModelDistance, DistanceRatio, CopyNumber, MajorChromosomeCount
        };
        public double GetQScorePredictor(QScorePredictor predictorId)
        {
            switch (predictorId)
            {
                case QScorePredictor.BinCount:
                    return (double)this.BinCount;

                case QScorePredictor.LogBinCount:
                    return Math.Log10(1 + this.BinCount);

                case QScorePredictor.BinCountAmpDistance:
                    return this.CopyNumber >= 15 ? Math.Log10(1 + this.BinCount) : 0.0;

                case QScorePredictor.BinMean:
                    if (this.Counts.Count == 0) return 0;
                    return this.Counts.Average();

                case QScorePredictor.BinCv:
                    if (this.Counts.Count == 0) return 0;
                    if (this.Counts.Average() == 0) return 0;
                    return Utilities.CoefficientOfVariation(this.Counts);

                case QScorePredictor.MafCount:
                    return Alleles.Frequencies.Count;

                case QScorePredictor.MafMean:
                    if (Alleles.Frequencies.Count == 0) return 0;
                    return Alleles.Frequencies.Average();

                case QScorePredictor.MafCv:
                    if (Alleles.Frequencies.Count == 0) return 0;
                    if (Alleles.Frequencies.Average() == 0) return 0;
                    return Utilities.CoefficientOfVariation(Alleles.Frequencies);

                case QScorePredictor.LogMafCv:
                    return Math.Log10(1 + GetQScorePredictor(QScorePredictor.MafCv));

                case QScorePredictor.ModelDistance:
                    return this.ModelDistance;

                case QScorePredictor.RunnerUpModelDistance:
                    return this.RunnerUpModelDistance;

                case QScorePredictor.DistanceRatio:
                    if (this.RunnerUpModelDistance == 0) return 0;
                    return this.ModelDistance / this.RunnerUpModelDistance;

                case QScorePredictor.CopyNumber:
                    return (double)this.CopyNumber;

                case QScorePredictor.MajorChromosomeCount:
                    // Force a double:
                    if (!this.MajorChromosomeCount.HasValue) return Math.Ceiling(this.CopyNumber / 2f);
                    return (double)this.MajorChromosomeCount;
            }
            return 0;
        }

        /// <summary>
        /// Set segment.Filter for each of our segments.
        /// </summary>
        public static void FilterSegments(int qualityFilterThreshold, List<CanvasSegment> segments)
        {
            string qualityFilter = $"q{qualityFilterThreshold}";
            foreach (var segment in segments)
            {
                string filter = null;
                if (segment.QScore < qualityFilterThreshold)
                {
                    filter = qualityFilter;
                }
                if (segment.End - segment.Begin < 10000)
                {
                    if (filter != null)
                        filter = filter + ";L10kb";
                    else
                        filter = "L10kb";
                }
                if (filter == null)
                    filter = "PASS";

                segment.Filter = filter;
            }
        }
    }
}
