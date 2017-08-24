using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using CanvasCommon;
using Illumina.Common;
using Isas.Framework.DataTypes;

namespace CanvasPedigreeCaller
{
    public class PedigreeMember
    {
        public enum Kinship
        {
            Other, Parent, Proband
        }
        
        public List<CanvasSegmentsSet> SegmentSets = new List<CanvasSegmentsSet>(); 
        

        public PedigreeMemberInfo PedigreeMemberInfo { get; set; }
        public SampleId Id;
        public string Name => Id.ToString();
        public Kinship Kin { get; set; }


        public double GetCoverage(int setPosition, int segmentPosition, SegmentsSet segmentsSet, int numberOfTrimmedBins)
        {
            return SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].MedianCount;
        }

        public List<Tuple<int, int>> GetAlleleCounts(int setPosition, int segmentPosition, SegmentsSet segmentsSet)
        {
            return SegmentSets[setPosition].GetSet(segmentsSet)[segmentPosition].Balleles.GetAlleleCounts();
        }
        public int GetPloidy(int haplotypeIndex, int segmentIndex, SegmentsSet segmentsSet)
        {
            return PedigreeMemberInfo.Ploidy?.GetReferenceCopyNumber(SegmentSets[haplotypeIndex].GetSet(segmentsSet)[segmentIndex]) ?? 2;
        }
        public PedigreeMemberInfo SetPedigreeMemberInfo(Segments segments, string ploidyBedPath, int numberOfTrimmedBins, Dictionary<string, List<Balleles>> allelesByChromosome,
            int maximumCopyNumber)
        {
            float meanMafCoverage = allelesByChromosome.SelectMany(x => x.Value).Select(y => y.MeanCoverage).Average();
            double variance = Utilities.Variance(segments.AllSegments.Select(x => x.TruncatedMedianCount(numberOfTrimmedBins)).ToList());
            double mafVariance = Utilities.Variance(segments.AllSegments.Where(x => x.Balleles.TotalCoverage.Count > 0)
                .Select(x => x.Balleles.TotalCoverage.Average()).ToList());
            double meanCoverage = segments.AllSegments.Select(x => x.TruncatedMedianCount(numberOfTrimmedBins)).Average();
            int maxCoverage = Convert.ToInt16(segments.AllSegments.Select(x => x.TruncatedMedianCount(numberOfTrimmedBins)).Max()) + 10;
            var ploidy = new PloidyInfo();
            if (!ploidyBedPath.IsNullOrEmpty() && File.Exists(ploidyBedPath))
                ploidy = PloidyInfo.LoadPloidyFromVcfFile(ploidyBedPath, Name);
            return new PedigreeMemberInfo(meanCoverage, meanMafCoverage, variance, mafVariance, maxCoverage, maximumCopyNumber, ploidy);
        }

    }
}