using System.Collections.Generic;
using System.IO;
using System.Linq;
using Canvas.CommandLineParsing;
using CanvasCommon;

namespace Canvas.SmallPedigree
{
    public class PedigreeSample
    {
        public SingleSampleCallset Sample { get; }
        public SampleType SampleType { get; }

        public PedigreeSample(SingleSampleCallset sample, SampleType sampleType)
        {
            Sample = sample;
            SampleType = sampleType;
        }
    }

    public class SmallPedigreeCallset
    {
        public AnalysisDetails AnalysisDetails { get; }
        public List<PedigreeSample> PedigreeSample;
        public bool HasPedigreeStructure { get; }
        public SmallPedigreeCallset(List<PedigreeSample> pedigreeSample, AnalysisDetails analysisDetails, bool hasPedigreeStructure)
        {
            PedigreeSample = pedigreeSample;
            AnalysisDetails = analysisDetails;
            HasPedigreeStructure = hasPedigreeStructure;
        }

        internal IEnumerable<string> NormalBinnedPath
        {
            get { return PedigreeSample.Select(sample => Path.Combine(sample.Sample.SampleOutputFolder.FullName, $"{sample.Sample.SampleName}.normal.binned")); }
        }

        internal IEnumerable<string> BinSizePath
        {
            get { return PedigreeSample.Select(sample => Path.Combine(sample.Sample.SampleOutputFolder.FullName, $"{sample.Sample.SampleName}.binsize")); }
        }

        internal IEnumerable<string> VfSummaryPath
        {
            get { return PedigreeSample.Select(sample => Path.Combine(sample.Sample.SampleOutputFolder.FullName, $"VFResults{sample.Sample.SampleName}.txt.gz")); }
        }
    }
}