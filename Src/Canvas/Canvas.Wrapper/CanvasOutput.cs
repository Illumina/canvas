using System;
using Illumina.Common.FileSystem;
using Isas.Framework.DataTypes;
using Isas.Framework.Utilities;

namespace Canvas.Wrapper
{
    public class CanvasOutput : ICanvasOutput
    {
        public Vcf CnvVcf { get; }
        public IFileLocation CoverageAndVariantFrequencies { get; }
        public IFileLocation VariantFrequencies { get; }
        public IFileLocation VariantFrequenciesBaf { get; }
        public IFileLocation Partitioned { get; }

        public CanvasOutput(
            Vcf cnvVcf,
            IFileLocation coverageAndVariantFrequencies,
            IFileLocation variantFrequencies = null,
            IFileLocation variantFrequenciesBaf = null,
            IFileLocation partitioned = null,
            IFileLocation binSize = null,
            IFileLocation normalBinned = null)
        {
            CnvVcf = cnvVcf;
            CoverageAndVariantFrequencies = coverageAndVariantFrequencies;
            VariantFrequencies = variantFrequencies;
            VariantFrequenciesBaf = variantFrequenciesBaf;
            Partitioned = partitioned;
        }

        public static CanvasOutput GetFromStub(IFileLocation stub)
        {
            Vcf cnvVcf = Vcf.GetVcfFromStub(stub);
            IFileLocation coverageAndVariantFrequencies = stub.AppendName(".CoverageAndVariantFrequency.txt");

            IFileLocation variantFrequencies = stub.AppendName(".VFResults.txt.gz");
            IFileLocation variantFrequenciesBaf = stub.AppendName(".VFResults.baf");
            IFileLocation partitioned = stub.AppendName(".partitioned");
            return new CanvasOutput(cnvVcf, coverageAndVariantFrequencies, variantFrequencies,
                variantFrequenciesBaf, partitioned);
        }

        public static string GetCoverageAndVariantFrequencyOutputPath(string outputVcfPath)
        {
            string coveragePath = outputVcfPath;
            if (outputVcfPath.EndsWith(".vcf.gz"))
                coveragePath = coveragePath.ReplaceFileNameExtension(".vcf.gz", "");
            else
                coveragePath.ReplaceFileNameExtension(".vcf", "");
            return coveragePath + ".CoverageAndVariantFrequency.txt";
        }

        public void Move(IFileLocation fileNameStub, Action<IFileLocation, IFileLocation> move)
        {
            CanvasOutput destination = GetFromStub(fileNameStub);
            CnvVcf.Move(destination.CnvVcf, move);
            move(CoverageAndVariantFrequencies, destination.CoverageAndVariantFrequencies);
            VariantFrequencies.MoveIfNotNull(destination.VariantFrequencies, move);
            VariantFrequenciesBaf.MoveIfNotNull(destination.VariantFrequenciesBaf, move);
            Partitioned.MoveIfNotNull(destination.Partitioned, move);
        }
    }
}