using Illumina.Common.FileSystem;
using Isas.Framework.DataTypes;

namespace Canvas.Wrapper.SmallPedigree
{
    public class IntermediateOutput
    {
        public IntermediateOutput(Vcf cnvVcf, IFileLocation coverageAndVariantFrequencies, IFileLocation variantFrequencies, IFileLocation variantFrequenciesBaf, 
            IFileLocation partitioned, IFileLocation coverageBigwig, BgzfFile bAlleleBedgraph, BgzfFile copyNumberBedgraph)
        {
            CnvVcf = cnvVcf;
#pragma warning disable CS0618 // Type or member is obsolete
            CoverageAndVariantFrequencies = coverageAndVariantFrequencies;
            VariantFrequencies = variantFrequencies;
            VariantFrequenciesBaf = variantFrequenciesBaf;
            Partitioned = partitioned;
#pragma warning restore CS0618 // Type or member is obsolete
            CoverageBigwig = coverageBigwig;
            BAlleleBedgraph = bAlleleBedgraph;
            CopyNumberBedgraph = copyNumberBedgraph;
        }

        public Vcf CnvVcf { get; }
        [System.Obsolete("Deprecated output - to be removed after NSv7 release")]
        public IFileLocation CoverageAndVariantFrequencies { get; }
        [System.Obsolete("Deprecated output - to be removed after NSv7 release")]
        public IFileLocation VariantFrequencies { get; }
        [System.Obsolete("Deprecated output - to be removed after NSv7 release")]
        public IFileLocation VariantFrequenciesBaf { get; }
        [System.Obsolete("Deprecated output - to be removed after NSv7 release")]
        public IFileLocation Partitioned { get; }

        public IFileLocation CoverageBigwig { get; }
        public BgzfFile BAlleleBedgraph { get; }
        public BgzfFile CopyNumberBedgraph { get; }

    }
}