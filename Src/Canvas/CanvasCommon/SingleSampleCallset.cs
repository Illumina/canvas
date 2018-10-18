using System.IO;
using Illumina.Common.FileSystem;
using Isas.Framework.DataTypes;

namespace CanvasCommon
{
    public class SingleSampleCallset
    {
        private readonly IDirectoryLocation _analysisOutputFolder;

        public SingleSampleCallset(Bam bam, string sampleName, IFileLocation normalVcfPath, bool isDbSnpVcf, IDirectoryLocation analysisOutputFolder, IFileLocation outputVcfPath)
        {
            _analysisOutputFolder = analysisOutputFolder;
            Bam = bam;
            SampleName = sampleName;
            NormalVcfPath = normalVcfPath;
            IsDbSnpVcf = isDbSnpVcf;
            OutputVcfPath = outputVcfPath;
        }

        public string SampleName { get; }
        public IFileLocation OutputVcfPath { get; }
        public Bam Bam { get; }
        public IFileLocation NormalVcfPath { get; }
        public bool IsDbSnpVcf { get; set; }
        public IDirectoryLocation SampleOutputFolder => GetSampleOutputFolder(_analysisOutputFolder, SampleName);
        public string BinSizePath => Path.Combine(SampleOutputFolder.FullName, $"{SampleName}.binsize");
        public string VfSummaryPath => GetVfSummaryPath(_analysisOutputFolder, SampleName).FullName;
        public string VfSummaryBafPath => GetVfSummaryBafPath(_analysisOutputFolder, SampleName).FullName;
        public BgzfFile BAlleleCoverageBedGraph => GetBAlleleBedGraph(_analysisOutputFolder, SampleName);
        public string NormalBinnedPath => Path.Combine(SampleOutputFolder.FullName, $"{SampleName}.normal.binned");
        public IFileLocation PartitionedPath => GetPartitionedPath(_analysisOutputFolder, SampleName);

        public static IDirectoryLocation GetSampleOutputFolder(IDirectoryLocation analysisOutputFolder, string sampleName)
        {
            return analysisOutputFolder.GetDirectoryLocation($"TempCNV_{sampleName}");
        }

        public static IFileLocation GetVfSummaryPath(IDirectoryLocation analysisOutputFolder, string sampleName)
        {
            var stub = GetSampleOutputFolder(analysisOutputFolder, sampleName).GetFileLocation($"VFResults{sampleName}");
            return GetVfSummaryPathExtension(stub);
        }

        public static IFileLocation GetVfSummaryBafPath(IDirectoryLocation analysisOutputFolder, string sampleName)
        {
            return GetVfSummaryPath(analysisOutputFolder, sampleName).AppendName(".baf");
        }

        public static IFileLocation GetVfSummaryBafPath(IFileLocation stub)
        {
            return stub.AppendName(".VFResults.baf");
        }

        private static IFileLocation GetVfSummaryPathExtension(IFileLocation stub)
        {
            return stub.AppendName(".txt.gz");
        }

        public static IFileLocation GetVfSummaryPath(IFileLocation stub)
        {
            return GetVfSummaryPathExtension(stub.AppendName(".VFResults"));
        }

        public static IFileLocation GetPartitionedPath(IDirectoryLocation analysisOutputFolder, string sampleName)
        {
            var stub = GetSampleOutputFolder(analysisOutputFolder, sampleName).GetFileLocation(sampleName);
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

        public static IFileLocation GetCoverageAndVariantFrequencyOutput(IDirectoryLocation analysisOutputFolder, string sampleName)
        {
            return GetCoverageAndVariantFrequencyOutput(GetSampleOutputFolder(analysisOutputFolder, sampleName).GetFileLocation("CNV"));
        }

        public static IFileLocation GetVcfOutput(IFileLocation stub)
        {
            return stub.AppendName(".vcf.gz");
        }

        public static IFileLocation GetVcfOutput(IDirectoryLocation analysisOutputFolder, string sampleName)
        {
            return GetVcfOutput(GetSampleOutputFolder(analysisOutputFolder, sampleName).GetFileLocation("CNV"));
        }

        public static IFileLocation GetCoverageAndVariantFrequencyOutputPath(string outputVcfPath)
        {
            string coveragePath = outputVcfPath;
            if (outputVcfPath.EndsWith(".vcf.gz"))
                coveragePath = coveragePath.ReplaceFileNameExtension(".vcf.gz", "");
            else
                coveragePath.ReplaceFileNameExtension(".vcf", "");
            return GetCoverageAndVariantFrequencyOutput(new FileLocation(coveragePath));
        }

        public static IFileLocation GetCoverageBigWig(IDirectoryLocation analysisOutputFolder, string sampleName)
        {
            return GetSampleOutputFolder(analysisOutputFolder, sampleName).GetFileLocation("coverage.bigWig");
        }

        public static IFileLocation GetCoverageBigWig(IFileLocation stub)
        {
            return stub.AppendName(".coverage.bigWig");
        }

        public static IFileLocation GetPartitionBedGraph(IDirectoryLocation analysisOutputFolder, string sampleName)
        {
            return GetSampleOutputFolder(analysisOutputFolder, sampleName).GetFileLocation("partition.bedgraph.gz");
        }

        public static BgzfFile GetPartitionBedGraph(IFileLocation stub)
        {
            return new BgzfFile(stub.AppendName(".partition.bedgraph.gz"));
        }

        public static BgzfFile GetCopyNumberBedGraph(IDirectoryLocation analysisOutputFolder, string sampleName)
        {
            return new BgzfFile(GetSampleOutputFolder(analysisOutputFolder, sampleName).GetFileLocation("copynumber.bedgraph.gz"));
        }

        public static BgzfFile GetCopyNumberBedGraph(IFileLocation stub)
        {
            return new BgzfFile(stub.AppendName(".copynumber.bedgraph.gz"));
        }

        public static BgzfFile GetBAlleleBedGraph(IDirectoryLocation analysisOutputFolder, string sampleName)
        {
            return new BgzfFile(GetSampleOutputFolder(analysisOutputFolder, sampleName).GetFileLocation("ballele.bedgraph.gz"));
        }

        public static BgzfFile GetBAlleleBedGraph(IFileLocation stub)
        {
            return new BgzfFile(stub.AppendName(".ballele.bedgraph.gz"));
        }

    }
}