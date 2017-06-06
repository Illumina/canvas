using Illumina.Common.FileSystem;
using Isas.Framework.DataTypes;

namespace Canvas.Wrapper.SmallPedigree
{
    public class IntermediateOutput
    {
        public IntermediateOutput(Vcf cnvVcf, IFileLocation coverageAndVariantFrequencies, IFileLocation variantFrequencies, IFileLocation variantFrequenciesBaf, IFileLocation partitioned)
        {
            CnvVcf = cnvVcf;
            CoverageAndVariantFrequencies = coverageAndVariantFrequencies;
            VariantFrequencies = variantFrequencies;
            VariantFrequenciesBaf = variantFrequenciesBaf;
            Partitioned = partitioned;
        }

        public Vcf CnvVcf { get; }
        public IFileLocation CoverageAndVariantFrequencies { get; }
        public IFileLocation VariantFrequencies { get; }
        public IFileLocation VariantFrequenciesBaf { get; }
        public IFileLocation Partitioned { get; }
    }
}