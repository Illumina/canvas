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
        public IDirectoryLocation OutputFolder => Callset.First().Callset.OutputFolder; 
        public IEnumerable<string> SampleNames { get { return Callset.Select(x => x.Callset.SampleName); } }
        public IEnumerable<Bam> BamPaths { get { return Callset.Select(x => x.Callset.Bam); } }
        public IFileLocation VcfPath => Callset.First().Callset.NormalVcfPath; 
        public IDirectoryLocation WholeGenomeFastaFolder { get { return Callset.Select(x => x.Callset.WholeGenomeFastaFolder).First(); } }
        public IFileLocation KmerFasta { get { return Callset.Select(x => x.Callset.KmerFasta).First(); } }
        public GenomeMetadata GenomeMetadata => Callset.First().Callset.GenomeMetadata; 
        public IFileLocation FilterBed => Callset.First().Callset.FilterBed; 
        public IFileLocation CommonCnvsBed { get; set; }
        public IFileLocation PloidyBed => Callset.First().Callset.PloidyVcf;
        public bool IsDbSnpVcf { get; set; } // NormalVcfPath points to a dbSNP VCF file

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