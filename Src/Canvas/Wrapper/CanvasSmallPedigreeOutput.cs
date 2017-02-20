using System;
using System.Collections.Generic;
using Illumina.Common.FileSystem;
using Isas.Framework.DataTypes;
using Isas.Framework.Utilities;

namespace Canvas.Wrapper
{
    public class CanvasSmallPedigreeTmpOutput
    {
        public CanvasSmallPedigreeTmpOutput(IFileLocation coverageAndVariantFrequencies, IFileLocation variantFrequencies, IFileLocation variantFrequenciesBaf, IFileLocation partitioned)
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

    public class CanvasSmallPedigreeOutput : ICanvasOutput
    {
        public List<CanvasSmallPedigreeTmpOutput> CanvasSmallPedigreeTmpOutput { get; set; }
        public Vcf CnvVcf { get; }

        public CanvasSmallPedigreeOutput(
            Vcf cnvVcf,
            List<IFileLocation> coverageAndVariantFrequencies,
            List<IFileLocation> variantFrequencies = null,
            List<IFileLocation> variantFrequenciesBaf = null,
            List<IFileLocation> partitioned = null)
        {
            CnvVcf = cnvVcf;
            var tmpOutputFiles = new List<CanvasSmallPedigreeTmpOutput>();
            for (int index = 0; index < coverageAndVariantFrequencies.Count; index++) {
                tmpOutputFiles.Add(new CanvasSmallPedigreeTmpOutput(coverageAndVariantFrequencies[index], variantFrequencies?[index], 
            variantFrequenciesBaf?[index], partitioned?[index]));
            }
            CanvasSmallPedigreeTmpOutput = tmpOutputFiles;
        }

        public static CanvasSmallPedigreeOutput GetFromStub(IFileLocation stub, bool includeIntermediateResults = false )
        {
            Vcf cnvVcf = Vcf.GetVcfFromStub(stub);
            // for now only copy multi-sample VCF
            return new CanvasSmallPedigreeOutput(cnvVcf, null);
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

        public void Move(IFileLocation fileNameStub, bool includeIntermediateResults, Action<IFileLocation, IFileLocation> move)
        {
            CanvasSmallPedigreeOutput destination = GetFromStub(fileNameStub, includeIntermediateResults);
            CnvVcf.Move(destination.CnvVcf, move);
        }
    }
}