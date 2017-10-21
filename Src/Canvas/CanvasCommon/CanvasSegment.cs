using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Isas.SequencingFiles;


namespace CanvasCommon
{
    public class Ballele
    {
        public int Position { get; }
        public int CountsA { get; }
        public int CountsB { get; }
        public float Frequency;
        public int TotalCoverage;

        private static float GetFrequency(int alleleACounts, int alleleBCounts)
        {
            return alleleBCounts / (float)GetTotalCoverage(alleleACounts, alleleBCounts);
        }

        public double GetMaxFrequency()
        {
            return Math.Max(CountsA, CountsB) / (double)GetTotalCoverage(CountsA, CountsB);
        }

        private static int GetTotalCoverage(int alleleACounts, int alleleBCounts)
        {
            return alleleACounts + alleleBCounts;
        }

        public Ballele(int position, int alleleACounts, int alleleBCounts)
        {
            Position = position;
            CountsA = alleleACounts;
            CountsB = alleleBCounts;
            Frequency = GetFrequency(alleleACounts, alleleBCounts);
            TotalCoverage = GetTotalCoverage(alleleACounts, alleleBCounts);
        }

        public Ballele(int position, float frequency, int totalCoverage, int alleleACounts, int alleleBCounts)
        {
            Position = position;
            Frequency = frequency;
            TotalCoverage = totalCoverage;
            CountsA = alleleACounts;
            CountsB = alleleBCounts;
        }
    }

    public class Balleles
    {

        private List<Ballele> Range;

        public Balleles()
        {
            Range = new List<Ballele>();
        }
        public Balleles(List<Ballele> alleles)
        {
            Range = alleles;
        }

        public void Add(Balleles alleles)
        {
            Range.AddRange(alleles.Range);
        }

        public void Add(Ballele allele)
        {
            Range.Add(allele);
        }

        public int Size()
        {
            return Range.Count;
        }

        public double? MeanMAF => Frequencies?.Average();
        public List<int> TotalCoverage => Range?.Select(allele => allele.TotalCoverage).ToList();

        public void PruneBalleles(double frequencyThreshold)
        {
            Range = Range?.Where(v => v.Frequency > frequencyThreshold).ToList();
        }

        public List<float> Frequencies => Range?.Select(allele => allele.Frequency).ToList();
        public List<double> MaxFrequencies => Range?.Select(allele => allele.GetMaxFrequency()).ToList();


        public List<Tuple<int, int>> GetAlleleCounts()
        {
            return Range.Select(allele => new Tuple<int, int>(allele.CountsA, allele.CountsB)).ToList();
        }

        public Tuple<int, int> MedianCounts(Balleles balleles)
        {
            var item1 = Utilities.Median(balleles.Range.Select(allele => Math.Max(allele.CountsA, allele.CountsB)).ToList());
            var item2 = Utilities.Median(balleles.Range.Select(allele => Math.Min(allele.CountsA, allele.CountsB)).ToList());
            return new Tuple<int, int>(item1, item2);
        }

        public Balleles GetBallelesSubrange(int start, int end, int defaultAlleleCountThreshold)
        {
            var array = Range.Where(x => x.Position >= start && x.Position <= end).ToList();
            return new Balleles(array);
        }

    }

    public class CoverageInfo
    {
        public Dictionary<string, uint[]> StartByChr = new Dictionary<string, uint[]>();
        public Dictionary<string, uint[]> EndByChr = new Dictionary<string, uint[]>();
        public Dictionary<string, double[]> CoverageByChr = new Dictionary<string, double[]>();
    }

    public enum SegmentsSet
    {
        SetA,
        SetB
    }

    /// <summary>
    /// A group of overlapping genomic intervals. 
    /// Used for keeping track of merged common and Canvas-called intervals when common CNV option is used.   
    /// </summary>
    public class OverlappingSegmentsRegion
    {
        public List<CanvasSegment> SetA { get; }
        public List<CanvasSegment> SetB { get; }

        protected SegmentsSet SelectedSet { get; set; }

        public List<CanvasSegment> GetSet()
        {
            return SelectedSet == SegmentsSet.SetA ? SetA : SetB;
        }

        public void SetSet(SegmentsSet set)
        {
            SelectedSet = set;
        }


        public OverlappingSegmentsRegion(List<CanvasSegment> setA, List<CanvasSegment> setB)
        {
            SetA = setA;
            SetB = setB;
            SelectedSet = SegmentsSet.SetA;
        }

        /// <summary>
        /// Used when there is just a single segment without any overlapping segments
        /// </summary>
        /// <param name="segment"></param>
        public OverlappingSegmentsRegion(CanvasSegment segment)
        {
            SetA = segment.ToEnumerable().ToList();
            SetB = null;
            SelectedSet = SegmentsSet.SetA;
        }
    }

    /// <summary>
    /// Contains information about a genomic interval. Has functions for computing copy numbers and their likelihoods.
    /// </summary>
    public partial class CanvasSegment
    {

        public CanvasSegment(string chr, int begin, int end, List<SampleGenomicBin> counts, Balleles balleles = null)
        {
            Chr = chr;
            GenomicBins = new List<SampleGenomicBin>();
            Begin = begin;
            End = end;
            GenomicBins = new List<SampleGenomicBin>(counts);
            CopyNumber = -1;
            SecondBestCopyNumber = -1;
            Balleles = balleles ?? new Balleles();
        }

        #region Members
        public List<SampleGenomicBin> GenomicBins;
        public int CopyNumber { get; set; }
        public int SecondBestCopyNumber { get; set; }
        public int? MajorChromosomeCount;
        public double QScore;
        public double? DqScore;
        public double? MajorChromosomeCountScore;
        public double ModelDistance;
        public double RunnerUpModelDistance;
        public bool CopyNumberSwapped;
        public bool IsCommonCnv;
        public bool IsHeterogeneous;
        private const int NumberVariantFrequencyBins = 100;
        private const int OverlapWindowThreshold = 500;
        public List<string> Filter = new List<string> {CanvasFilter.Pass};
        public Tuple<int, int> StartConfidenceInterval; // if not null, this is a confidence interval around Start, reported in the CIPOS tag
        public Tuple<int, int> EndConfidenceInterval; // if not null, this is a confidence interval around End, reported in the CIEND tag
        public Balleles Balleles;
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
        public List<float> Counts
        {
            get
            {
                return this.GenomicBins.Select(x => x.Count).ToList();
            }
        }
        #endregion
        public int Length => End - Begin;
        public List<SampleGenomicBin> GetSampleGenomicBinSubrange(int start, int end)
        {
            return this.GenomicBins.Where(x => x.Start >= start && x.Stop <= end).ToList();
        }

        /// <summary>
        /// Mean of the segment's counts.
        /// </summary>
        public double MeanCount
        {
            get
            {
                double sum = Counts.Sum();
                return sum / this.BinCount;
            }
        }

        /// <summary>
        /// Median of the segment's counts.
        /// </summary>
        public double MedianCount
        {
            get
            {
                var sorted = new SortedList<double>(Counts.Select(Convert.ToDouble));
                return sorted.Median();
            }
        }

        /// <summary>
        /// removes flanking bins before median estimation
        /// </summary>
        public double TruncatedMedianCount(int bins2Remove)
        {
            int start = bins2Remove;
            int end = Counts.Count - bins2Remove;
            int length = end - start;
            var sorted = length > 5 ?
            new SortedList<double>(Counts.Skip(start).Take(length).Select(Convert.ToDouble)) :
            new SortedList<double>(Counts.Select(Convert.ToDouble));
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
            this.GenomicBins.AddRange(s.GenomicBins);
            if (this.Balleles != null && s.Balleles != null) Balleles.Add(s.Balleles);
        }

        public int SizeOveralp(CanvasSegment segment)
        {
            if (segment.Begin > this.Begin && segment.End < this.End)
                return this.Length - segment.Length;
            if (segment.Begin > this.Begin && segment.Begin < this.End && this.End <= segment.End)
                return this.End - segment.Begin;
            if (segment.Begin < this.Begin && segment.End > this.Begin && this.End > segment.End)
                return segment.End - this.Begin;
            return 0;
        }

        /// <summary>
        /// Handle different overlap scenarios between Canvas segment and common segment
        /// </summary>
        /// <param name="canvasSegments"></param>
        /// <param name="commonSegments"></param>
        /// <param name="defaultAlleleCountThreshold"></param>
        /// <param name="canvasSegmentsIndex"></param>
        /// <param name="commonSegmentsIndex"></param>
        /// <returns></returns>
        public static OverlappingSegmentsRegion SplitCanvasSegments(List<CanvasSegment> canvasSegments, List<CanvasSegment> commonSegments, int defaultAlleleCountThreshold,
            ref int canvasSegmentsIndex, ref int commonSegmentsIndex)
        {
            var haplotypebSegments = new List<CanvasSegment>();
            var haplotypeaSegments = new List<CanvasSegment>();

            // scenario: common segment within Canvas segment
            // canvasSegment:   ----------------------------------
            // commonSegment:         -----------------
            if (commonSegments[commonSegmentsIndex].Begin > canvasSegments[canvasSegmentsIndex].Begin &&
                commonSegments[commonSegmentsIndex].End < canvasSegments[canvasSegmentsIndex].End)
            {
                int begin = canvasSegments[canvasSegmentsIndex].Begin;
                int end = commonSegments[commonSegmentsIndex].Begin;
                var countsSubRange = canvasSegments[canvasSegmentsIndex].GetSampleGenomicBinSubrange(begin, end);
                var allelesSubRange = canvasSegments[canvasSegmentsIndex].Balleles.GetBallelesSubrange(begin, end, defaultAlleleCountThreshold);
                if (!countsSubRange.Empty())
                    haplotypebSegments.Add(new CanvasSegment(commonSegments[commonSegmentsIndex].Chr, begin, end, countsSubRange, allelesSubRange));
                haplotypebSegments.Add(commonSegments[commonSegmentsIndex]);

                // subscenario: Canvas segment spans more than one common segment
                // canvasSegment:   ------------------------------------------------
                // commonSegment:            ------------     -------------------
                if (commonSegments.Count > commonSegmentsIndex + 1 && commonSegments[commonSegmentsIndex + 1].Begin < canvasSegments[canvasSegmentsIndex].End)
                {
                    commonSegmentsIndex++;

                    while (commonSegments.Count > commonSegmentsIndex && commonSegments[commonSegmentsIndex].Begin < canvasSegments[canvasSegmentsIndex].End)
                    {
                        haplotypebSegments.Add(commonSegments[commonSegmentsIndex]);
                        commonSegmentsIndex++;
                    }
                    haplotypeaSegments.Add(canvasSegments[canvasSegmentsIndex]);
                    canvasSegmentsIndex++;
                    return new OverlappingSegmentsRegion(haplotypeaSegments, haplotypebSegments);
                }

                begin = commonSegments[commonSegmentsIndex].End;
                end = canvasSegments[canvasSegmentsIndex].End;
                countsSubRange = canvasSegments[canvasSegmentsIndex].GetSampleGenomicBinSubrange(begin, end);
                allelesSubRange = canvasSegments[canvasSegmentsIndex].Balleles.GetBallelesSubrange(begin, end, defaultAlleleCountThreshold);
                if (!countsSubRange.Empty())
                    haplotypebSegments.Add(new CanvasSegment(commonSegments[commonSegmentsIndex].Chr, begin, end, countsSubRange, allelesSubRange));

                haplotypeaSegments.Add(canvasSegments[canvasSegmentsIndex]);
                canvasSegmentsIndex++;
                commonSegmentsIndex++;
                return new OverlappingSegmentsRegion(haplotypeaSegments, haplotypebSegments);
            }

            // scenario: Canvas segment part-overlaps common segment and comes first
            // canvasSegment:   --------------------
            // commonSegment:            ------------------
            if (commonSegments[commonSegmentsIndex].Begin > canvasSegments[canvasSegmentsIndex].Begin && commonSegments[commonSegmentsIndex].Begin < canvasSegments[canvasSegmentsIndex].End &&
                canvasSegments[canvasSegmentsIndex].End <= commonSegments[commonSegmentsIndex].End)
            {
                haplotypeaSegments.Add(canvasSegments[canvasSegmentsIndex]);
                int begin = canvasSegments[canvasSegmentsIndex].Begin;
                int end = commonSegments[commonSegmentsIndex].Begin;
                var countsSubRange = canvasSegments[canvasSegmentsIndex].GetSampleGenomicBinSubrange(begin, end);
                var allelesSubRange = canvasSegments[canvasSegmentsIndex].Balleles.GetBallelesSubrange(begin, end, defaultAlleleCountThreshold);
                if (!countsSubRange.Empty())
                {
                    haplotypebSegments.Add(new CanvasSegment(commonSegments[commonSegmentsIndex].Chr, begin, end, countsSubRange, allelesSubRange));
                }

                // scenario: Canvas segment part-overlaps common segment,  comes first and both segments have the same end coords
                // canvasSegment:   ---------------------------
                // commonSegment:            ------------------
                if (canvasSegments[canvasSegmentsIndex].End == commonSegments[commonSegmentsIndex].End)
                {
                    haplotypebSegments.Add(commonSegments[commonSegmentsIndex]);
                    canvasSegmentsIndex++;
                    commonSegmentsIndex++;
                    return new OverlappingSegmentsRegion(haplotypeaSegments, haplotypebSegments);
                }

                // scenario: common segment spans more than one Canvas segment
                // canvasSegment:   -------------------      --------
                // commonSegment:            ---------------------------------
                if (canvasSegments.Count > canvasSegmentsIndex + 1 && commonSegments[commonSegmentsIndex].End > canvasSegments[canvasSegmentsIndex + 1].End)
                {
                    canvasSegmentsIndex++;

                    while (canvasSegments.Count > canvasSegmentsIndex && commonSegments[commonSegmentsIndex].End > canvasSegments[canvasSegmentsIndex].End)
                    {
                        haplotypeaSegments.Add(canvasSegments[canvasSegmentsIndex]);
                        canvasSegmentsIndex++;
                    }
                    haplotypebSegments.Add(commonSegments[commonSegmentsIndex]);
                    commonSegmentsIndex++;
                    return new OverlappingSegmentsRegion(haplotypeaSegments, haplotypebSegments);
                }

                haplotypebSegments.Add(commonSegments[commonSegmentsIndex]);
                canvasSegmentsIndex++;
                begin = canvasSegments[canvasSegmentsIndex].Begin;
                end = commonSegments[commonSegmentsIndex].End;
                countsSubRange = canvasSegments[canvasSegmentsIndex].GetSampleGenomicBinSubrange(begin, end);
                allelesSubRange = canvasSegments[canvasSegmentsIndex].Balleles.GetBallelesSubrange(begin, end, defaultAlleleCountThreshold);
                if (!countsSubRange.Empty())
                    haplotypeaSegments.Add(new CanvasSegment(commonSegments[commonSegmentsIndex].Chr, begin, end, countsSubRange, allelesSubRange));
                canvasSegments[commonSegmentsIndex].Begin = commonSegments[commonSegmentsIndex].End + 1;
                return new OverlappingSegmentsRegion(haplotypeaSegments, haplotypebSegments);
            }

            // scenario: common segment part-overlaps Canvas segment and comes first
            // canvasSegment:           --------------------
            // commonSegment:   ------------------
            if (commonSegments[commonSegmentsIndex].Begin <= canvasSegments[canvasSegmentsIndex].Begin && commonSegments[commonSegmentsIndex].End >
                canvasSegments[canvasSegmentsIndex].Begin && canvasSegments[canvasSegmentsIndex].End > commonSegments[commonSegmentsIndex].End)
            {
                haplotypebSegments.Add(commonSegments[commonSegmentsIndex]);

                int begin = commonSegments[commonSegmentsIndex].End;
                int end = canvasSegments[canvasSegmentsIndex].End;
                var countsSubRange = canvasSegments[canvasSegmentsIndex].GetSampleGenomicBinSubrange(begin, end);
                var allelesSubRange = canvasSegments[canvasSegmentsIndex].Balleles.GetBallelesSubrange(begin, end, defaultAlleleCountThreshold);
                if (!countsSubRange.Empty())
                    haplotypebSegments.Add(new CanvasSegment(commonSegments[commonSegmentsIndex].Chr, begin, end, countsSubRange, allelesSubRange));

                haplotypeaSegments.Add(canvasSegments[canvasSegmentsIndex]);

                canvasSegmentsIndex++;
                commonSegmentsIndex++;
                return new OverlappingSegmentsRegion(haplotypeaSegments, haplotypebSegments);
            }

            // default: for now do now handle other conditions
            canvasSegmentsIndex++;
            return new OverlappingSegmentsRegion(new List<CanvasSegment> { canvasSegments[canvasSegmentsIndex] }, null);
        }

        public int BinCount => Counts.Count;

        /// <summary>
        /// Compute the median count from a list of segments.
        /// </summary>
        /// <param name="segments">List of segments.</param>
        /// <returns>Median of counts.</returns>
        public static double ExpectedCount(List<CanvasSegment> segments)
        {
            var counts = new List<double>();

            // Get all of the counts contained within all autosomal segments.
            foreach (CanvasSegment segment in segments)
            {
                if (GenomeMetadata.SequenceMetadata.IsAutosome(segment.Chr))
                {
                    counts.AddRange(segment.Counts.Select(count => (double)count));
                }
            }
            double mu = Utilities.Median(counts);
            return mu;
        }

        /// <summary>
        /// Apply quality scores.
        /// </summary>
        public static void AssignQualityScores(IReadOnlyList<CanvasSegment> segments, QScoreMethod qscoreMethod, QualityScoreParameters qscoreParameters)
        {
            foreach (CanvasSegment segment in segments)
            {
                segment.QScore = segment.ComputeQScore(qscoreMethod, qscoreParameters);
            }
        }

        public static ConcurrentDictionary<string, List<CanvasSegment>> GetSegmentsByChromosome(IReadOnlyList<CanvasSegment> segments)
        {
            var segmentsByChromosome = new ConcurrentDictionary<string, List<CanvasSegment>>();
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

        public static Dictionary<string, List<CanvasSegment>> GetCommonCnvSegmentsByChromosome(List<CanvasSegment> segments)
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
        private static int GetMinimumBinsForCoveragePlotPoint(IReadOnlyList<CanvasSegment> segments, int pointLength)
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
        public static void WriteCoveragePlotData(IReadOnlyList<CanvasSegment> segments, double? normalDiploidCoverage, PloidyInfo referencePloidy,
            IFileLocation file, string referenceFolder)
        {
            int pointLength = 100000;
            int minimumBinsToPlot = GetMinimumBinsForCoveragePlotPoint(segments, pointLength);

            var segmentsByChromosome = GetSegmentsByChromosome(segments);
            GenomeMetadata genome = new GenomeMetadata();
            genome.Deserialize(new FileLocation(Path.Combine(referenceFolder, "GenomeSize.xml")));


            List<float> counts = new List<float>();
            List<float> MAF = new List<float>();
            List<float> VF = new List<float>();
            file.Directory.Create();
            using (FileStream stream = new FileStream(file.FullName, FileMode.Create, FileAccess.Write))
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
                            string key = String.Format("{0} {1}", segment.CopyNumber, segment.MajorChromosomeCount);
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
                            if (Int32.Parse(bits[0]) != majorCopyNumber) continue;
                            long count = CopyNumberAndChromCount[key];
                            if (count < bestCount) continue;
                            bestCount = count;
                            majorChromosomeCount = Int32.Parse(bits[1]);
                        }

                        // Note allele frequency and coverage info, for all overlap segments that match (more or less)
                        // the most common copy number:
                        foreach (CanvasSegment segment in overlapSegments)
                        {
                            if ((majorCopyNumber == 2 && segment.CopyNumber != 2) ||
                                (majorCopyNumber < 2 && segment.CopyNumber >= 2) ||
                                (majorCopyNumber > 2 && segment.CopyNumber <= 2))
                                continue;

                            // Add counts to the overall list:
                            int firstIndex = 0;
                            if (pointStartPos > segment.Begin)
                            {
                                firstIndex = (int)((float)segment.GenomicBins.Count * (pointStartPos - segment.Begin) / segment.Length);
                            }
                            int lastIndex = segment.Counts.Count;
                            if (pointEndPos < segment.End)
                            {
                                lastIndex = (int)((float)segment.GenomicBins.Count * (pointEndPos - segment.Begin) / segment.Length);
                            }
                            for (int index = firstIndex; index < lastIndex; index++) counts.Add(segment.GenomicBins[index].Count);

                            // Add MAF to the overall list:
                            if (segment.Balleles != null)
                            {
                                firstIndex = 0;
                                if (pointStartPos > segment.Begin)
                                {
                                    firstIndex = (int)((float)segment.Balleles.Size() * (pointStartPos - segment.Begin) / segment.Length);
                                }
                                lastIndex = segment.Balleles.Size();
                                if (pointEndPos < segment.End)
                                {
                                    lastIndex = (int)((float)segment.Balleles.Size() * (pointEndPos - segment.Begin) / segment.Length);
                                }
                                VF.AddRange(segment.Balleles.Frequencies.Skip(firstIndex).Take(lastIndex - firstIndex));
                                MAF.AddRange(segment.Balleles.MaxFrequencies.Select(Convert.ToSingle).Skip(firstIndex).Take(lastIndex - firstIndex));
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
        private static bool IsForbiddenInterval(string chr, int start, int end,
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
        /// Merge common CNV regions and Canvas segment creating CanvasSegmentsHaplotype blocks of CanvasSegments.
        /// CanvasSegmentsHaplotype will have two haplotypes where common CNV regions and Canvas segment overlap
        /// The method calculates various overlap scenarious and outputs list of CanvasSegmentsHaplotype
        /// </summary>
        /// <param name="canvasSegments"></param>
        /// <param name="commonCnvSegments"></param>
        /// <param name="chr"></param>
        /// <param name="defaultAlleleCountThreshold"></param>
        /// <returns></returns>
        public static List<OverlappingSegmentsRegion> MergeCommonCnvSegments(List<CanvasSegment> canvasSegments,
            List<CanvasSegment> commonCnvSegments, int defaultAlleleCountThreshold)
        {
            const int segmentOverlapThreshold = 10;
            var mergedSegments = new List<OverlappingSegmentsRegion>(canvasSegments.Count + commonCnvSegments.Count * 3);
            var sortedCanvasSegments = canvasSegments.OrderBy(o => o.Begin).ToList();
            var sortedCommonCnvSegments = commonCnvSegments.OrderBy(o => o.Begin).ToList();
            var canvasSegmentsIndex = 0;
            var commonSegmentsIndex = 0;
            if (sortedCanvasSegments[canvasSegmentsIndex].End <= sortedCommonCnvSegments[commonSegmentsIndex].Begin)
            {
                mergedSegments.Add(new OverlappingSegmentsRegion(new List<CanvasSegment> { sortedCanvasSegments[0] }, null));
                canvasSegmentsIndex++;
            }

            // iterate over two CanvasSegment lists and merge using various scenarious 
            while (canvasSegmentsIndex < sortedCanvasSegments.Count && commonSegmentsIndex < sortedCommonCnvSegments.Count)
            {

                // skip small common CNV variants 
                if (sortedCommonCnvSegments[commonSegmentsIndex].Length < OverlapWindowThreshold * 2)
                {
                    commonSegmentsIndex++;
                    continue;
                }
                // Canvas segment comes before common segment and does not overlap it
                if (sortedCanvasSegments[canvasSegmentsIndex].End <= sortedCommonCnvSegments[commonSegmentsIndex].Begin)
                {
                    mergedSegments.Add(new OverlappingSegmentsRegion(new List<CanvasSegment> { sortedCanvasSegments[canvasSegmentsIndex] }, null));
                    canvasSegmentsIndex++;
                    continue;
                }
                // common segment comes before Canvas segment and does not overlap it
                if (sortedCanvasSegments[canvasSegmentsIndex].Begin >= sortedCommonCnvSegments[commonSegmentsIndex].End)
                {
                    mergedSegments.Add(new OverlappingSegmentsRegion(null, new List<CanvasSegment> { sortedCommonCnvSegments[commonSegmentsIndex] }));
                    commonSegmentsIndex++;
                    continue;
                }
                // Canvas segment and common segment have the same coordinates
                if (sortedCanvasSegments[canvasSegmentsIndex].Begin == sortedCommonCnvSegments[commonSegmentsIndex].Begin &&
                    sortedCanvasSegments[canvasSegmentsIndex].End == sortedCommonCnvSegments[commonSegmentsIndex].End)
                {
                    mergedSegments.Add(new OverlappingSegmentsRegion(null, new List<CanvasSegment> { sortedCommonCnvSegments[commonSegmentsIndex] }));
                    canvasSegmentsIndex++;
                    commonSegmentsIndex++;
                    continue;
                }
                // Canvas segment and common segment have coordinates within segmentation margin of error
                if (Math.Abs(sortedCanvasSegments[canvasSegmentsIndex].Begin - sortedCommonCnvSegments[commonSegmentsIndex].Begin) < OverlapWindowThreshold &&
                    Math.Abs(sortedCanvasSegments[canvasSegmentsIndex].End - sortedCommonCnvSegments[commonSegmentsIndex].End) < OverlapWindowThreshold &&
                    sortedCommonCnvSegments[commonSegmentsIndex].Length > OverlapWindowThreshold * 4)
                {
                    mergedSegments.Add(new OverlappingSegmentsRegion(null, new List<CanvasSegment> { sortedCommonCnvSegments[commonSegmentsIndex] }));
                    canvasSegmentsIndex++;
                    commonSegmentsIndex++;
                    continue;
                }
                // common segment and Canvas segment overlap
                if (sortedCanvasSegments[canvasSegmentsIndex].SizeOveralp(sortedCommonCnvSegments[commonSegmentsIndex]) > segmentOverlapThreshold)
                {
                    var newSegmentsHaplotype = SplitCanvasSegments(sortedCanvasSegments, sortedCommonCnvSegments, defaultAlleleCountThreshold,
                        ref canvasSegmentsIndex, ref commonSegmentsIndex);
                    mergedSegments.Add(newSegmentsHaplotype);
                }
                else
                {
                    mergedSegments.Add(new OverlappingSegmentsRegion(new List<CanvasSegment> { sortedCanvasSegments[canvasSegmentsIndex] }, null));
                    canvasSegmentsIndex++;
                    commonSegmentsIndex++;
                }
            }

            if (canvasSegmentsIndex < sortedCanvasSegments.Count)
                mergedSegments.AddRange(sortedCanvasSegments.Skip(canvasSegmentsIndex).Select(segment => new OverlappingSegmentsRegion(new List<CanvasSegment> { segment }, null)));

            else if (commonSegmentsIndex < sortedCommonCnvSegments.Count)
                mergedSegments.AddRange(sortedCommonCnvSegments.Skip(commonSegmentsIndex).Select(segment => new OverlappingSegmentsRegion(null, new List<CanvasSegment> { segment })));

            return mergedSegments;
        }

        /// <summary>
        /// Iterates through a list of segments and merges those which have the same copy number call.
        /// Also, for segments smaller than MinimumCallSize, assimilate them into the neighbor with the best 
        /// quality score.  Two consecutive segments are considered neighbors if they're on the same chromosome
        /// and the space between them doesn't overlap with any excluded intervals.
        /// </summary>
        public static List<CanvasSegment> MergeSegmentsUsingExcludedIntervals(List<CanvasSegment> segments, int MinimumCallSize,
            Dictionary<string, List<SampleGenomicBin>> excludedIntervals)
        {
            // Assimilate short segments into the *best* available neighbor:
            var mergedSegments = new List<CanvasSegment>();
            if (!segments.Any()) return mergedSegments;

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
            return mergedSegments;
        }

        /// <summary>
        /// Iterates through a list of segments and merges those which have the same copy number call.
        /// For multisample workflow a 2D list of regions x samples is provided to test for identity of CN
        /// calls across all samples.
        /// Also, for segments smaller than MinimumCallSize, assimilate them into the neighbor with the best 
        /// quality score.  Two consecutive segments are considered neighbors if they're on the same chromosome
        /// and the space between them is not too large.
        /// </summary>
        public static List<CanvasSegment> MergeSegments(List<CanvasSegment> segments, int minimumCallSize = 0, int maximumMergeSpan = 10000,
            List<List<int>> copyNumbers = null, List<double> qscores = null)
        {
            // Assimilate short segments into the *best* available neighbor:
            var mergedSegments = new List<CanvasSegment>();
            if (!segments.Any()) return mergedSegments;

            var newCopyNumbers = new List<List<int>>();
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
            return mergedSegments;
        }

        /// <summary>
        /// Set segment.Filter for each of our segments.
        /// </summary>
        public static void SetFilterForSegments(int qualityFilterThreshold, List<CanvasSegment> segments, int segmantSizeCutoff)
        {
            string qualityFilter = $"q{qualityFilterThreshold}";
            string sizeFilter = "L" + CanvasFilter.FormatCnvSizeWithSuffix(segmantSizeCutoff);
            foreach (var segment in segments)
            {
                if (segment.Filter.Count > 0) throw new Exception($"Filter has already been set: {string.Join(";", segment.Filter)}");
                if (segment.QScore < qualityFilterThreshold) segment.Filter.Add(qualityFilter);
                if (segment.End - segment.Begin < segmantSizeCutoff) segment.Filter.Add(sizeFilter);
                if (segment.Filter.Count == 0) segment.Filter.Add(CanvasFilter.Pass); // "PASS" is exclusive with any other keyword
            }
        }

        /// <summary>
        /// Remap GenomicBin from genome coordiantes into CanvasBin coordiantes
        /// TODO: Stoping using this method. We don't want to use indexes and we certainly don't want to create bed intervals to represent ranges of indexes 
        /// </summary>
        public static List<BedInterval> RemapGenomicToBinCoordinates(List<BedEntry> commonRegions, IReadOnlyList<SampleGenomicBin> sampleGenomicBins)
        {
            var searchStartIndex = 0;
            var commonRegionsRemapped = new List<BedInterval>();
            var intervals = sampleGenomicBins.Select(bin => bin.GenomicBin.Interval).ToList();
            foreach (var commonRegion in commonRegions)
            {
                int genomicBinStartIndex = intervals.FindIndex(searchStartIndex, interval => interval.Start <= commonRegion.Start && interval.End > commonRegion.Start);
                int genomicBinEndIndex = intervals.FindIndex(searchStartIndex, interval => interval.Start <= commonRegion.End && interval.End > commonRegion.End);
                if (genomicBinStartIndex == -1 || genomicBinEndIndex == -1) continue;
                commonRegionsRemapped.Add(new BedInterval(genomicBinStartIndex, genomicBinEndIndex));
                searchStartIndex = genomicBinEndIndex;
            }
            return commonRegionsRemapped;
        }


        /// <summary>
        /// Assume that the rows are sorted by the start position and ascending order
        /// </summary>
        public static CoverageInfo ReadBEDInput(string inputbed, string forbiddenIntervalBedPath = null)
        {
            const int idxChr = 0, idxStart = 1, idxEnd = 2, idxScore = 3;
            var binFilter = new GenomicBinFilter(forbiddenIntervalBedPath);
            var coverageInfo = new CoverageInfo();
            try
            {
                var startByChr = new Dictionary<string, List<uint>>();
                var endByChr = new Dictionary<string, List<uint>>();
                var scoreByChr = new Dictionary<string, List<double>>();
                using (var reader = new GzipReader(inputbed))
                {
                    string line;
                    while ((line = reader.ReadLine()) != null)
                    {
                        var tokens = line.Split('\t');
                        string chrom = tokens[idxChr].Trim();
                        uint start = Convert.ToUInt32(tokens[idxStart].Trim());
                        uint end = Convert.ToUInt32(tokens[idxEnd].Trim());
                        if (binFilter.SkipBin(chrom, start, end))
                            continue;
                        if (!startByChr.ContainsKey(chrom))
                        {
                            startByChr.Add(chrom, new List<uint>());
                            endByChr.Add(chrom, new List<uint>());
                            scoreByChr.Add(chrom, new List<double>());
                        }
                        startByChr[chrom].Add(start);
                        endByChr[chrom].Add(end);
                        scoreByChr[chrom].Add(Convert.ToDouble(tokens[idxScore].Trim()));
                    }
                    foreach (string chr in startByChr.Keys)
                    {
                        coverageInfo.StartByChr[chr] = startByChr[chr].ToArray();
                        coverageInfo.EndByChr[chr] = endByChr[chr].ToArray();
                        coverageInfo.CoverageByChr[chr] = scoreByChr[chr].ToArray();
                    }
                }
            }
            catch (Exception e)
            {
                Console.Error.WriteLine("File {0} could not be read:", inputbed);
                Console.Error.WriteLine(e.Message);
                Environment.Exit(1);
            }
            return coverageInfo;
        }

        public static List<CanvasSegment> CreateSegmentsFromCommonCnvs(IReadOnlyList<SampleGenomicBin> sampleGenomicBins, List<BedInterval> commonSegmentIntervals,
            List<Balleles> allelesForCommonSegments)
        {
            var segments = new List<CanvasSegment>();
            if (commonSegmentIntervals.Last().End > sampleGenomicBins.Count)
                throw new IndexOutOfRangeException("Coverage bin index exceeds chromosome size (in Canvas bins)");
            foreach (var seg in commonSegmentIntervals.Zip(allelesForCommonSegments, (s, a) => (interval: s, alleles: a)))
            {
                int length = seg.interval.End - seg.interval.Start;
                var segment = new CanvasSegment(sampleGenomicBins[seg.interval.Start].GenomicBin.Chromosome,
                    sampleGenomicBins[seg.interval.Start].Start,
                    sampleGenomicBins[seg.interval.End].Stop,
                    sampleGenomicBins.Skip(seg.interval.Start).Take(length).ToList(),
                    seg.alleles)
                {
                    IsCommonCnv = true
                };
                segments.Add(segment);
            }
            return segments;
        }

    }
}
