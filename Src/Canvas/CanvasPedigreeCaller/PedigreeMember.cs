using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using CanvasCommon;
using Illumina.Common;

namespace CanvasPedigreeCaller
{
    public class PedigreeMemberInfo
    {
        public PedigreeMemberInfo(double meanCoverage, double meanMafCoverage, double variance, double mafVariance, int maxCoverage, int maximumCopyNumber, PloidyInfo ploidy)
        {
            MeanCoverage = meanCoverage;
            MeanMafCoverage = meanMafCoverage;
            Variance = variance;
            MafVariance = mafVariance;
            MaxCoverage = maxCoverage;
            Ploidy = ploidy;
            CnModel = new CopyNumberModel(maximumCopyNumber, this);
        }

        public double MeanCoverage { get; internal set; }
        public double MeanMafCoverage { get; internal set; }
        public double Variance { get; internal set; }
        public double MafVariance { get; internal set; }
        public int MaxCoverage { get; internal set; }
        public PloidyInfo Ploidy { get; internal set; }
        public CopyNumberModel CnModel { get; set; }

    }

    public class PedigreeMember
    {
        public enum Kinship
        {
            Other, Parent, Proband
        }

        public Segments SegmentsByChromosome;
        public List<CanvasSegment> Segments = new List<CanvasSegment>();
        public List<SegmentHaplotypes> SegmentSets = new List<SegmentHaplotypes>();
        public PedigreeMemberInfo PedigreeMemberInfo { get; set; }
        public string Name { get; set; }
        public List<string> Parents { get; set; }
        public List<string> Offspring { get; set; }
        public Kinship Kin { get; set; }


        public double GetCoverage(int setPosition, int segmentPosition, Haplotype haplotype, int numberOfTrimmedBins)
        {
            return SegmentSets[setPosition].GetSet(haplotype)[segmentPosition].MedianCount;
        }

        public List<Tuple<int, int>> GetAlleleCounts(int setPosition, int segmentPosition, Haplotype haplotype)
        {
            return SegmentSets[setPosition].GetSet(haplotype)[segmentPosition].Balleles.GetAlleleCounts();
        }
        public int GetPloidy(int haplotypeIndex, int segmentIndex, Haplotype haplotype)
        {
            return PedigreeMemberInfo.Ploidy?.GetReferenceCopyNumber(SegmentSets[haplotypeIndex].GetSet(haplotype)[segmentIndex]) ?? 2;
        }
        public PedigreeMemberInfo SetPedigreeMemberInfo(string ploidyBedPath, int numberOfTrimmedBins, Dictionary<string, List<Balleles>> allelesByChromosome, 
            int maximumCopyNumber)
        {
            float meanMafCoverage = allelesByChromosome.SelectMany(x => x.Value).Select(y => y.MeanCoverage).Average();
            double variance = Utilities.Variance(Segments.Select(x => x.TruncatedMedianCount(numberOfTrimmedBins)).ToList());
            double mafVariance = Utilities.Variance(Segments.Where(x => x.Balleles.TotalCoverage.Count > 0)
                .Select(x => x.Balleles.TotalCoverage.Average()).ToList());
            double meanCoverage = Segments.Select(x => x.TruncatedMedianCount(numberOfTrimmedBins)).Average();
            int maxCoverage = Convert.ToInt16(Segments.Select(x => x.TruncatedMedianCount(numberOfTrimmedBins)).Max()) + 10;
            var ploidy = new PloidyInfo();
            if (!ploidyBedPath.IsNullOrEmpty() && File.Exists(ploidyBedPath))
                ploidy = PloidyInfo.LoadPloidyFromVcfFile(ploidyBedPath, Name);
            return new PedigreeMemberInfo(meanCoverage, meanMafCoverage, variance, mafVariance, maxCoverage, maximumCopyNumber, ploidy);
        }

    }
}