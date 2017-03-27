using Illumina.Common.FileSystem;

namespace Canvas.Wrapper.SmallPedigree
{
    public class IntermediateOutput
    {
        public IntermediateOutput(IFileLocation coverageAndVariantFrequencies, IFileLocation variantFrequencies, IFileLocation variantFrequenciesBaf, IFileLocation partitioned)
        {
            CoverageAndVariantFrequencies = coverageAndVariantFrequencies;
            VariantFrequencies = variantFrequencies;
            VariantFrequenciesBaf = variantFrequenciesBaf;
            Partitioned = partitioned;
        }

        public IFileLocation CoverageAndVariantFrequencies { get; }
        public IFileLocation VariantFrequencies { get; }
        public IFileLocation VariantFrequenciesBaf { get; }
        public IFileLocation Partitioned { get; }
    }
}