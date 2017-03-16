using Isas.Framework.DataTypes;

namespace Canvas.Wrapper.SmallPedigree
{
    public class CanvasSmallPedigreeOutput
    {
        public SampleSet<IntermediateOutput> IntermediateOutputs { get; set; }
        public Vcf CnvVcf { get; }

        public CanvasSmallPedigreeOutput(
            Vcf cnvVcf,
            SampleSet<IntermediateOutput> intermediateOutputs)
        {
            CnvVcf = cnvVcf;
            IntermediateOutputs = intermediateOutputs;
        }
    }
}