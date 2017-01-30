namespace Illumina.SecondaryAnalysis.VariantCalling.StructuralVariants.Canvas
{
    public class CanvasTumorNormalWgsInput : ICanvasCheckpointInput
    {
        public Bam TumorBam { get; }
        public Bam NormalBam { get; }
        public Vcf NormalVcf { get; } // set to the starling VCF path 
        public Vcf SomaticVcf { get; } // set to the strelka VCF path
        public GenomeMetadata GenomeMetadata { get; }

        public CanvasTumorNormalWgsInput(Bam tumorBam, Bam normalBam, Vcf normalVcf, Vcf somaticVcf, GenomeMetadata genomeMetadata)
        {
            TumorBam = tumorBam;
            NormalBam = normalBam;
            NormalVcf = normalVcf;
            SomaticVcf = somaticVcf;
            GenomeMetadata = genomeMetadata;
        }
    }
}