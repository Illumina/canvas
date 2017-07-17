using System;
using System.Collections.Generic;
using System.Linq;
using CanvasCommon;

namespace CanvasPedigreeCaller
{
    public class PedigreeMember
    {
        public enum Kinship
        {
            Parent, Offspring, Proband
        }
        public List<CanvasSegment> Segments = new List<CanvasSegment>();
        public List<CanvasSegmentsSet> SegmentSets = new List<CanvasSegmentsSet>();
        public double MeanCoverage { get; set; }
        public int MaxCoverage { get; set; }
        public double MeanMafCoverage { get; set; }
        public string Name { get; set; }
        public List<string> Parents { get; set; }
        public List<string> Offspring { get; set; }
        public PloidyInfo Ploidy { get; set; }
        public double Variance { get; internal set; }
        public double MafVariance { get; internal set; }
        public CopyNumberModel CnModel { get; set; }
        public Kinship Kin { get; set; }

        public double GetCoverage(int setPosition, int segmentPosition, SegmentsSet segmentsSet, int numberOfTrimmedBins)
        {
            return SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].MedianCount;
        }
        public Tuple<int,int> GetMedianAlleleCounts(int segmentIndex)
        {
            int allele1 = Math.Min(Segments[segmentIndex].Alleles.MedianCounts.Item1, MaxCoverage - 1);
            int allele2 = Math.Min(Segments[segmentIndex].Alleles.MedianCounts.Item2, MaxCoverage - 1);
            return new Tuple<int, int>(allele1, allele2);
        }
        public List<Tuple<int, int>> GetAlleleCounts(int setPosition, int segmentPosition, SegmentsSet segmentsSet)
        {
            return SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].Alleles.Balleles.Select(x=>x.Counts).ToList();
        }

        public int GetPloidy(int haplotypeIndex, int segmentIndex, SegmentsSet segmentsSet)
        {
            return Ploidy?.GetReferenceCopyNumber(SegmentSets[haplotypeIndex].GetSet(segmentsSet)[segmentIndex]) ?? 2;
        }
    }
}