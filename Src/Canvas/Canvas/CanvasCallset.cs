using System.Collections.Generic;
using System.IO;
using System.Linq;
using Illumina.Common.FileSystem;
using Isas.Framework.DataTypes;
using Isas.Manifests.NexteraManifest;
using Isas.SequencingFiles;

namespace Canvas
{

    public class AnalysisDetails
    {
        public IDirectoryLocation OutputFolder { get; }
        public IDirectoryLocation WholeGenomeFastaFolder { get; }
        public IFileLocation KmerFasta { get;  }
        public GenomeMetadata GenomeMetadata { get; }
        public IFileLocation FilterBed { get; }
        public IFileLocation PloidyVcf { get; }
        public IFileLocation CommonCnvsBed { get; }

        public AnalysisDetails(IDirectoryLocation outputFolder, IDirectoryLocation wholeGenomeFastaFolder,
           IFileLocation kmerFasta, IFileLocation filterBed, IFileLocation ploidyVcf, IFileLocation commonCnvsBed)        
        {
            WholeGenomeFastaFolder = wholeGenomeFastaFolder;
            OutputFolder = outputFolder;
            KmerFasta = kmerFasta;
            FilterBed = filterBed;
            PloidyVcf = ploidyVcf;
            CommonCnvsBed = commonCnvsBed;
            var genomeSizeXml = WholeGenomeFastaFolder.GetFileLocation("GenomeSize.xml");
            GenomeMetadata = new GenomeMetadata();
            GenomeMetadata.Deserialize(genomeSizeXml.FullName);
        }
        internal string TempFolder => Path.Combine(OutputFolder.FullName, "TempCNV");
    }


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
        internal string TempFolder => Path.Combine(OutputFolder.FullName, $"TempCNV_{SampleName}");
        internal string BinSizePath => Path.Combine(TempFolder, $"{SampleName}.binsize");
        internal string VfSummaryPath => Path.Combine(TempFolder, $"VFResults{SampleName}.txt.gz");
        internal string VfSummaryBafPath => VfSummaryPath + ".baf";
        internal string NormalBinnedPath => Path.Combine(TempFolder, $"{SampleName}.normal.binned");
    }

    public class CanvasCallset
    {
        public readonly SingleSampleCallset SingleSampleCallset;
        public IEnumerable<Bam> NormalBamPaths { get; }
        public IFileLocation SomaticVcfPath { get; } // set to the strelka VCF path
        public AnalysisDetails AnalysisDetails { get; set; }
        public NexteraManifest Manifest { get; }


        public CanvasCallset(
            IFileLocation bam,
            string sampleName,
            IFileLocation normalVcfPath,
            bool isDbSnpVcf,
            IEnumerable<IFileLocation> normalBamPaths,
            NexteraManifest manifest,
            IFileLocation somaticVcfPath,
            IFileLocation outputVcfPath,
            AnalysisDetails analysisDetails) 

        {
            SingleSampleCallset = new SingleSampleCallset(new Bam(bam), sampleName, normalVcfPath, isDbSnpVcf, analysisDetails.OutputFolder, outputVcfPath);
            Manifest = manifest;
            SomaticVcfPath = somaticVcfPath;
            AnalysisDetails = analysisDetails;
            NormalBamPaths = normalBamPaths.Select(file => new Bam(file));
        }

        public CanvasCallset(
            SingleSampleCallset singleSampleCallset,
            AnalysisDetails analysisDetails,
            IEnumerable<IFileLocation> normalBamPaths,
            NexteraManifest manifest,
            IFileLocation somaticVcfPath)

        {
            SingleSampleCallset = singleSampleCallset;
            Manifest = manifest;
            if (somaticVcfPath != null)
                SomaticVcfPath = somaticVcfPath;
            AnalysisDetails = analysisDetails;
            if (normalBamPaths != null)
                NormalBamPaths = normalBamPaths.Select(file => new Bam(file));
        }
        public bool IsEnrichment => Manifest != null;
        internal string TempManifestPath => Path.Combine(SingleSampleCallset.TempFolder, "manifest.txt");
    }
}