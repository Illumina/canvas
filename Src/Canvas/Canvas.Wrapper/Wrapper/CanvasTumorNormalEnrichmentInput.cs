using Illumina.SecondaryAnalysis.VariantCalling;
using Isas.Framework.DataTypes;
using Isas.Manifests.NexteraManifest;
using Isas.SequencingFiles;

namespace Canvas.Wrapper
{
    public class CanvasTumorNormalEnrichmentInput : ICanvasEnrichmentInput
    {
        public Bam TumorBam { get; }
        public Bam NormalBam { get; }
        public Vcf NormalVcf { get; } // set to the Starling VCF path (if tumor normal, the normal vcf path) 
        public Vcf SomaticVcf { get; } // set to the strelka VCF path
        public GenomeMetadata GenomeMetadata { get; }
        public NexteraManifest NexteraManifest { get; }
        public SexPloidyInfo SexPloidy { get; }

        public CanvasTumorNormalEnrichmentInput(
            Bam tumorBam,
            Bam normalBam,
            Vcf normalVcf,
            Vcf somaticVcf,
            GenomeMetadata genomeMetadata,
            NexteraManifest nexteraManifest, SexPloidyInfo sexPloidy)
        {
            TumorBam = tumorBam;
            NormalBam = normalBam;
            NormalVcf = normalVcf;
            SomaticVcf = somaticVcf;
            GenomeMetadata = genomeMetadata;
            NexteraManifest = nexteraManifest;
            SexPloidy = sexPloidy;
        }
    }
}