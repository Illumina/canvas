using Illumina.SecondaryAnalysis.VariantCalling;
using Isas.Framework.DataTypes;
using Isas.SequencingFiles;

namespace Canvas.Wrapper
{
    public class CanvasResequencingInput : ICanvasCheckpointInput
    {
        public Bam Bam { get; }
        public GenomeMetadata GenomeMetadata { get; }
        public SexPloidyInfo SexPloidy { get; }
        public Vcf Vcf { get; }

        public CanvasResequencingInput(Bam bam, Vcf vcf, GenomeMetadata genomeMetadata, SexPloidyInfo sexPloidy)
        {
            Bam = bam;
            GenomeMetadata = genomeMetadata;
            SexPloidy = sexPloidy;
            Vcf = vcf;
        }
    }
}