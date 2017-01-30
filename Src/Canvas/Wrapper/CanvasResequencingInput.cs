namespace Illumina.SecondaryAnalysis.VariantCalling.StructuralVariants.Canvas
{
    public class CanvasResequencingInput : ICanvasCheckpointInput
    {
        public Bam Bam { get; }
        public GenomeMetadata GenomeMetadata { get; }
        public Vcf Vcf { get; }

        public CanvasResequencingInput(Bam bam, Vcf vcf, GenomeMetadata genomeMetadata)
        {
            Bam = bam;
            GenomeMetadata = genomeMetadata;
            Vcf = vcf;
        }
    }
}