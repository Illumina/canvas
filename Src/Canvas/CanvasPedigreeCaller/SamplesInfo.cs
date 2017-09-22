using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using CanvasCommon;
using Isas.Framework.DataTypes;
using Illumina.Common;

namespace CanvasPedigreeCaller
{
    public class SamplesInfo
    {

        private SamplesInfo(double meanCoverage, double meanMafCoverage, double variance, 
            double mafVariance, int maxCoverage, PloidyInfo ploidy, SampleId id)
        {
            MeanCoverage = meanCoverage;
            MeanMafCoverage = meanMafCoverage;
            Variance = variance;
            MafVariance = mafVariance;
            MaxCoverage = maxCoverage;
            Ploidy = ploidy;
            SampleId = id;
        }

        public double MeanCoverage { get; }
        public SampleId SampleId { get; }
        public double MeanMafCoverage { get; }
        public double Variance { get; }
        public double MafVariance { get; }
        public int MaxCoverage { get; }
        public PloidyInfo Ploidy { get; }

        public int GetPloidy(CanvasSegment segment)
        {
            return Ploidy?.GetReferenceCopyNumber(segment) ?? 2;
        }
        public static SamplesInfo GetPedigreeMemberInfo(Segments segments, string ploidyBedPath, int numberOfTrimmedBins, SampleId id)
        {
            double meanMafCoverage = segments.AllSegments.SelectMany(x => x.Balleles.TotalCoverage).Average();
            double variance = Utilities.Variance(segments.AllSegments.Select(x => x.TruncatedMedianCount(numberOfTrimmedBins)).ToList());
            double mafVariance = Utilities.Variance(segments.AllSegments.Where(x => x.Balleles.TotalCoverage.Count > 0)
                .Select(x => x.Balleles.TotalCoverage.Average()).ToList());
            double meanCoverage = segments.AllSegments.Select(x => x.TruncatedMedianCount(numberOfTrimmedBins)).Average();
            int maxCoverage = Convert.ToInt16(segments.AllSegments.Select(x => x.TruncatedMedianCount(numberOfTrimmedBins)).Max()) + 10;
            var ploidy = new PloidyInfo();
            if (!ploidyBedPath.IsNullOrEmpty() && File.Exists(ploidyBedPath))
                ploidy = PloidyInfo.LoadPloidyFromVcfFile(ploidyBedPath, id.ToString());
            return new SamplesInfo(meanCoverage, meanMafCoverage, variance, mafVariance, maxCoverage, ploidy, id);
        }
    }
}
