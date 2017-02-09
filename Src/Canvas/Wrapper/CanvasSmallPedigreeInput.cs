using Isas.Framework.DataTypes;
using Isas.SequencingFiles;

namespace Canvas.Wrapper
{
    public class CanvasSmallPedigreeInput : ICanvasCheckpointInput
    {
        public Bam Bam { get; }
        public GenomeMetadata GenomeMetadata { get; }
        public Vcf Vcf { get; }

        public CanvasSmallPedigreeInput(Bam bam, Vcf vcf, GenomeMetadata genomeMetadata)
        {
            Bam = bam;
            GenomeMetadata = genomeMetadata;
            Vcf = vcf;
        }
    }
}