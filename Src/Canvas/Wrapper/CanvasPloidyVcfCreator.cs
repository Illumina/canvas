using System.Linq;
using Illumina.Common.FileSystem;
using Illumina.SecondaryAnalysis.VariantCalling;
using Isas.Framework.DataTypes;
using Isas.Framework.Logging;
using Isas.Framework.WorkManagement;
using Isas.SequencingFiles;

namespace Canvas.Wrapper
{
    public class CanvasPloidyVcfCreator
    {
        public CanvasPloidyVcfCreator(ILogger logger, IWorkManager workManager, PloidyCorrector ploidyFixer)
        {
            _logger = logger;
            _workManager = workManager;
            _ploidyFixer = ploidyFixer;
        }

        private readonly IWorkManager _workManager;
        private readonly ILogger _logger;
        private readonly PloidyCorrector _ploidyFixer;

        /// <summary>
        /// Write out the ploidy bed file if ploidy information is available from the vcf header
        /// </summary>
        public IFileLocation CreatePloidyVcf(SampleSet<SexPloidyInfo> ploidyInfos, GenomeMetadata genomeMetadata, IDirectoryLocation sampleSandbox)
        {
            IFileLocation ploidyVcf = sampleSandbox.GetFileLocation("ploidy.vcf.gz");
            _ploidyFixer.WritePloidyVcfFile()
            string fastaPath = genomeMetadata.Sequences.First().FastaPath;
            if (_ploidyFixer.GeneratePloidyBedFileFromVcf(
                genomeMetadata,
                fastaPath,
                vcf.VcfFile.FullName,
                ploidyBed.FullName, sampleSandbox.FullName, _logger, _workManager))
            {
                return ploidyBed;
            }
            _logger.Warn($"Sex chromosome ploidy not found in {vcf.VcfFile} header. No ploidy will be provided to Canvas.");
            return null;
        }
    }
}