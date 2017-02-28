using System.IO;
using Illumina.Common.FileSystem;
using Isas.Framework.DataTypes;

namespace CanvasCommon
{
    public class SingleSampleCallset
    {
        public SingleSampleCallset(Bam bam, string sampleName, IFileLocation normalVcfPath, bool isDbSnpVcf, IDirectoryLocation outputFolder, IFileLocation outputVcfPath)
        {
            Bam = bam;
            SampleName = sampleName;
            NormalVcfPath = normalVcfPath;
            IsDbSnpVcf = isDbSnpVcf;
            OutputFolder = outputFolder;
            OutputVcfPath = outputVcfPath;
        }

        public string SampleName { get; }
        public IFileLocation OutputVcfPath { get; }
        public Bam Bam { get; }
        public IFileLocation NormalVcfPath { get; }
        public bool IsDbSnpVcf { get; set; }
        public IDirectoryLocation OutputFolder { get; }
        private IDirectoryLocation TempDirectory => GetSampleTempFolder(OutputFolder, SampleName);
        public string TempFolder => TempDirectory.FullName;
        public string BinSizePath => Path.Combine(TempFolder, $"{SampleName}.binsize");
        public string VfSummaryPath => GetVfSummaryPath(TempDirectory, SampleName).FullName;
        public string VfSummaryBafPath => GetVfSummaryBafPath(TempDirectory, SampleName).FullName;
        public string NormalBinnedPath => Path.Combine(TempFolder, $"{SampleName}.normal.binned");
        public IFileLocation PartitionedPath => GetPartitionedPath(OutputFolder, SampleName);

        public static IDirectoryLocation GetSampleTempFolder(IDirectoryLocation outputFolder, string sampleName)
        {
            return outputFolder.GetDirectoryLocation($"TempCNV_{sampleName}");
        }

        public static IFileLocation GetVfSummaryPath(IDirectoryLocation outputFolder, string sampleName)
        {
            var stub = GetSampleTempFolder(outputFolder, sampleName).GetFileLocation($"VFResults{sampleName}");
            return GetVfSummaryPathExtension(stub);
        }

        public static IFileLocation GetVfSummaryBafPath(IDirectoryLocation outputFolder, string sampleName)
        {
            return GetVfSummaryPath(outputFolder, sampleName).AppendName(".baf");
        }

        public static IFileLocation GetVfSummaryBafPath(IFileLocation stub)
        {
            return GetVfSummaryPath(stub).AppendName(".baf");
        }

        private static IFileLocation GetVfSummaryPathExtension(IFileLocation stub)
        {
            return stub.AppendName(".txt.gz");
        }

        public static IFileLocation GetVfSummaryPath(IFileLocation stub)
        {
            return GetVfSummaryPathExtension(stub.AppendName(".VFResults"));
        }

        public static IFileLocation GetPartitionedPath(IDirectoryLocation outputFolder, string sampleName)
        {
            var stub = GetSampleTempFolder(outputFolder, sampleName).GetFileLocation(sampleName);
            return GetPartitionedPath(stub);
        }

        public static IFileLocation GetPartitionedPath(IFileLocation stub)
        {
            return stub.AppendName(".partitioned");
        }

        public static IFileLocation GetCoverageAndVariantFrequencyOutput(IFileLocation stub)
        {
            return stub.AppendName(".CoverageAndVariantFrequency.txt");
        }

        public static IFileLocation GetCoverageAndVariantFrequencyOutput(IDirectoryLocation output, string sampleName)
        {
            return GetCoverageAndVariantFrequencyOutput(SingleSampleCallset.GetSampleTempFolder(output, sampleName).GetFileLocation("CNV"));
        }

        public static string GetCoverageAndVariantFrequencyOutputPath(string outputVcfPath)
        {
            string coveragePath = outputVcfPath;
            if (outputVcfPath.EndsWith(".vcf.gz"))
                coveragePath = coveragePath.ReplaceFileNameExtension(".vcf.gz", "");
            else
                coveragePath.ReplaceFileNameExtension(".vcf", "");
            return GetCoverageAndVariantFrequencyOutput(new FileLocation(coveragePath)).FullName;
        }
    }
}