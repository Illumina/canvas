using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Illumina.SecondaryAnalysis;
using Isas.Shared;
using SequencingFiles;

namespace Canvas
{
    public class CanvasCallset
    {
        public IDirectoryLocation WholeGenomeFastaFolder { get; }
        public IDirectoryLocation OutputFolder { get; }
        public IFileLocation KmerFasta { get; }
        public string Id => SampleName; // unique ID for this callset
        public string SampleName { get; }
        public Bam Bam { get; }
        public IEnumerable<Bam> NormalBamPaths { get; }
        public IFileLocation NormalVcfPath { get; }// set to the Starling VCF path (if tumor normal, the normal vcf path) 
        public bool IsDbSnpVcf { get; set; }// NormalVcfPath points to a dbSNP VCF file
        public IFileLocation SomaticVcfPath { get; } // set to the strelka VCF path
        public IFileLocation OutputVcfPath { get; }
        public NexteraManifest Manifest { get; }
        public GenomeMetadata GenomeMetadata { get; }
        public IFileLocation FilterBed { get; }
        public IFileLocation PloidyBed { get; }

        public CanvasCallset(
            IFileLocation bam,
            string sampleName,
            IDirectoryLocation wholeGenomeFastaFolder,
            IDirectoryLocation outputFolder,
            IFileLocation kmerFasta,
            IFileLocation filterBed,
            IFileLocation ploidyBed,
            IFileLocation normalVcfPath,
            bool isDbSnpVcf,
            IEnumerable<IFileLocation> normalBamPaths,
            NexteraManifest manifest,
            IFileLocation somaticVcfPath,
            IFileLocation outputVcfPath)
        {
            Bam = new Bam(bam);
            SampleName = sampleName;
            WholeGenomeFastaFolder = wholeGenomeFastaFolder;
            OutputFolder = outputFolder;
            KmerFasta = kmerFasta;
            FilterBed = filterBed;
            PloidyBed = ploidyBed;
            NormalVcfPath = normalVcfPath;
            IsDbSnpVcf = isDbSnpVcf;
            Manifest = manifest;
            SomaticVcfPath = somaticVcfPath;
            OutputVcfPath = outputVcfPath;
            NormalBamPaths = normalBamPaths.Select(file => new Bam(file));

            var genomeSizeXml = WholeGenomeFastaFolder.GetFileLocation("GenomeSize.xml");
            GenomeMetadata = new GenomeMetadata();
            GenomeMetadata.Deserialize(genomeSizeXml.FullName);
        }

        public bool IsEnrichment => Manifest != null;

        internal string TempFolder
        {
            get { return Path.Combine(OutputFolder.FullName, String.Format("TempCNV_{0}", Id)); }
        }

        internal string NormalBinnedPath
        {
            get { return Path.Combine(TempFolder, String.Format("{0}.normal.binned", Id)); }
        }

        internal string BinSizePath
        {
            get { return Path.Combine(TempFolder, String.Format("{0}.binsize", Id)); }
        }

        internal string VfSummaryPath
        {
            get { return Path.Combine(TempFolder, String.Format("VFResults{0}.txt.gz", Id)); }
        }

        internal string TempManifestPath
        {
            get { return Path.Combine(TempFolder, "manifest.txt"); }
        }
    }
}