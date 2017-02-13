using Illumina.SecondaryAnalysis.VariantCalling;
using Isas.Framework.DataTypes;
using Isas.SequencingFiles;

namespace Canvas.Wrapper
{
    public class CanvasSmallPedigreeInput : ICanvasCheckpointInput
    {
        public Bam Bam { get; }
        public GenomeMetadata GenomeMetadata { get; }
        public Vcf Vcf { get; }
        public SampleSet<SexPloidyInfo> PloidyInfos { get; set; }
        public string PedigreeName { get; set; }

        public CanvasSmallPedigreeInput(Bam bam, Vcf vcf, GenomeMetadata genomeMetadata)
        {
            Bam = bam;
            GenomeMetadata = genomeMetadata;
            Vcf = vcf;
        }
    }
}