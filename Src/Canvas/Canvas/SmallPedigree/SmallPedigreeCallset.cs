using System.Collections.Generic;
using System.IO;
using System.Linq;
using Canvas.CommandLineParsing;
using Isas.SequencingFiles;
using Isas.Shared.DataTypes;
using Isas.Shared.Utilities.FileSystem;

namespace Canvas.SmallPedigree
{
    public class SingleSampleCallset
    {
        public CanvasCallset Callset { get; }
        public SampleType SampleType { get; }

        public SingleSampleCallset(CanvasCallset callset, SampleType sampleType)
        {
            Callset = callset;
            SampleType = sampleType;
        }
    }
    public class SmallPedigreeCallset
    {
        public IDirectoryLocation OutputFolder { get { return Callset.Select(x => x.Callset.OutputFolder).First(); } }
        public IEnumerable<string> SampleNames { get { return Callset.Select(x => x.Callset.SampleName); } }
        public IEnumerable<Bam> BamPaths { get { return Callset.Select(x => x.Callset.Bam); } }
        public IEnumerable<IFileLocation> NormalVcfPaths { get; } // set to the Starling VCF path (if tumor normal, the normal vcf path) 
        public IDirectoryLocation WholeGenomeFastaFolder { get; set; }
        public IFileLocation KmerFasta { get { return Callset.Select(x => x.Callset.KmerFasta).First(); } }
        public GenomeMetadata GenomeMetadata { get { return Callset.Select(x => x.Callset.GenomeMetadata).First(); } }
        public IFileLocation FilterBed { get { return Callset.Select(x => x.Callset.FilterBed).First(); } }
        public IFileLocation CommonCnvsBed { get; set; }
        public IFileLocation PloidyBed { get { return Callset.Select(x => x.Callset.PloidyBed).First(); } }
        public bool IsDbSnpVcf { get; set; } // NormalVcfPath points to a dbSNP VCF file
        public IFileLocation OutputVcfPath { get; }
        public NexteraManifest Manifest { get; }

        public List<SingleSampleCallset> Callset;
        public SmallPedigreeCallset(List<SingleSampleCallset> callset, IFileLocation commonCnvsBed)
        {
            Callset = callset;
            CommonCnvsBed = commonCnvsBed;
        }

        internal string TempFolder
        {
            get { return Path.Combine(OutputFolder.FullName, "TempCNV"); }
        }

        internal IEnumerable<string> NormalBinnedPath
        {
            get { return SampleNames.Select(sampleName => Path.Combine(TempFolder, $"{sampleName}.normal.binned")); }
        }

        internal IEnumerable<string> BinSizePath
        {
            get { return SampleNames.Select(sampleName => Path.Combine(TempFolder, $"{sampleName}.binsize")); }
        }

        internal IEnumerable<string> VfSummaryPath
        {
            get { return SampleNames.Select(sampleName => Path.Combine(TempFolder, $"VFResults{sampleName}.txt.gz")); }
        }
    }
}