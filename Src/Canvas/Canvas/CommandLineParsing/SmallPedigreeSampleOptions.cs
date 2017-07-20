using Illumina.Common.FileSystem;

namespace Canvas.CommandLineParsing
{
    public class SmallPedigreeSampleOptions
    {
        public string SampleName { get; }
        public SampleType SampleType { get; }
        public IFileLocation Bam { get; }

        public SmallPedigreeSampleOptions(string sampleName, SampleType sampleType, IFileLocation bam)
        {
            SampleName = sampleName;
            SampleType = sampleType;
            Bam = bam;
        }
    }
}