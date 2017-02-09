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
        public Vcf CreatePloidyVcf(SampleSet<SexPloidyInfo> ploidyInfos, GenomeMetadata genomeMetadata, IDirectoryLocation sampleSandbox)
        {
            var ploidyVcf = new Vcf(sampleSandbox.GetFileLocation("ploidy.vcf.gz"));
            _ploidyFixer.WritePloidyVcfFile(ploidyVcf, ploidyInfos, genomeMetadata);
            return ploidyVcf;
        }
    }
}