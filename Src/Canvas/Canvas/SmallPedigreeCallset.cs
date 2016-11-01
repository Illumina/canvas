using System.Collections.Generic;
using System.IO;
using System.Linq;
using Isas.SequencingFiles;
using Isas.Shared.DataTypes;
using Isas.Shared.Utilities.FileSystem;

namespace Canvas
{
    public class SmallPedigreeCallset
    {
        public IDirectoryLocation OutputFolder { get { return _singleSampleCallset.Select(x => x.OutputFolder).First(); } }
        public IEnumerable<string> SampleNames { get { return _singleSampleCallset.Select(x => x.SampleName); } }
        public IEnumerable<Bam> BamPaths { get { return _singleSampleCallset.Select(x => x.Bam); } }
        public IEnumerable<IFileLocation> NormalVcfPaths { get; } // set to the Starling VCF path (if tumor normal, the normal vcf path) 
        public IDirectoryLocation WholeGenomeFastaFolder { get; set; }
        public IFileLocation KmerFasta { get { return _singleSampleCallset.Select(x => x.KmerFasta).First(); } }
        public GenomeMetadata GenomeMetadata { get { return _singleSampleCallset.Select(x => x.GenomeMetadata).First(); } }
        public IFileLocation FilterBed { get { return _singleSampleCallset.Select(x => x.FilterBed).First(); } }
        public IFileLocation CommonCnvsBed { get; set; }
        public IFileLocation PedigreeInfo { get; set; }
        public IFileLocation PloidyBed { get { return _singleSampleCallset.Select(x => x.PloidyBed).First(); } }
        public bool IsDbSnpVcf { get; set; } // NormalVcfPath points to a dbSNP VCF file
        public IFileLocation OutputVcfPath { get; }
        public NexteraManifest Manifest { get; }

        private readonly List<CanvasCallset> _singleSampleCallset;
        public SmallPedigreeCallset(List<CanvasCallset> callset, IFileLocation commonCnvsBed, IFileLocation pedigreeInfo)
        {
            _singleSampleCallset = callset;
            CommonCnvsBed = commonCnvsBed;
            PedigreeInfo = pedigreeInfo;
        }

        public List<CanvasCallset> Callset { get { return _singleSampleCallset; } }

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