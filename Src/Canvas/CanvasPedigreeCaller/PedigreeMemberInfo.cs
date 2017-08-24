using System.Collections.Generic;
using System.Linq;
using CanvasCommon;
using Isas.Framework.DataTypes;

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

        public double MeanCoverage { get; }
        public double MeanMafCoverage { get; }
        public double Variance { get; }
        public double MafVariance { get; }
        public int MaxCoverage { get; }
        public PloidyInfo Ploidy { get; }
        public CopyNumberModel CnModel { get; }
    }
}