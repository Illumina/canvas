using Canvas.CommandLineParsing;
using CanvasCommon;
using Illumina.SecondaryAnalysis.VariantCalling;
using Isas.Framework.DataTypes;
using Isas.SequencingFiles;

namespace Canvas.Wrapper.SmallPedigree
{
    public class CanvasSmallPedigreeInput : ICanvasCheckpointInput
    {
        public GenomeMetadata GenomeMetadata { get; }
        public Vcf Vcf { get; }
        public SampleSet<CanvasPedigreeSample> Samples { get; }

        public CanvasSmallPedigreeInput(GenomeMetadata genomeMetadata, SampleSet<CanvasPedigreeSample> samples, Vcf vcf)
        {
            GenomeMetadata = genomeMetadata;
            Samples = samples;
            Vcf = vcf;
        }
    }

    public class CanvasPedigreeSample
    {
        public CanvasPedigreeSample(Bam bam, SampleType sampleType, SexPloidyInfo ploidyInfo)
        {
            Bam = bam;
            SampleType = sampleType;
            PloidyInfo = ploidyInfo;
        }

        public Bam Bam { get; }
        public SampleType SampleType { get; }
        public SexPloidyInfo PloidyInfo { get; }
    }
}