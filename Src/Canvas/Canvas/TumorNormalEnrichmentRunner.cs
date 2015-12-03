using System.Linq;
using Canvas.CommandLineParsing;
using CanvasCommon;
using Illumina.SecondaryAnalysis;
using Illumina.SecondaryAnalysis.Workflow;
using Isas.Shared;
using SequencingFiles;

namespace Canvas
{
    public class TumorNormalEnrichmentRunner : IModeRunner
    {
        private readonly TumorNormalOptions _tumorNormalOptions;
        private readonly IFileLocation _normalBam;
        private readonly IFileLocation _manifest;
        public CommonOptions CommonOptions { get; }

        public TumorNormalEnrichmentRunner(CommonOptions commonOptions, TumorNormalOptions tumorNormalOptions, IFileLocation normalBam, IFileLocation manifest)
        {
            _tumorNormalOptions = tumorNormalOptions;
            _normalBam = normalBam;
            _manifest = manifest;
            CommonOptions = commonOptions;
        }

        public void Run(ILogger logger, ICheckpointRunner checkpointRunner, IWorkManager workManager)
        {
            CanvasRunner runner = new CanvasRunner(logger, workManager, checkpointRunner, false, CanvasCoverageMode.GCContentWeighted, 300, CommonOptions.CustomParams);
            var callset = GetCallset(logger);
            runner.CallSample(callset);
        }

        private CanvasCallset GetCallset(ILogger logger)
        {
            IFileLocation outputVcfPath = CommonOptions.OutputDirectory.GetFileLocation("CNV.vcf.gz");
            var manifest = new NexteraManifest(_manifest.FullName, null, logger.Error);
            CanvasCallset callSet = new CanvasCallset(
                    _tumorNormalOptions.TumorBam,
                    CommonOptions.SampleName,
                    CommonOptions.WholeGenomeFasta,
                    CommonOptions.OutputDirectory,
                    CommonOptions.KmerFasta,
                    CommonOptions.FilterBed,
                    CommonOptions.PloidyBed,
                    CommonOptions.BAlleleSites,
                    CommonOptions.IsDbSnpVcf,
                    new[] { _normalBam },
                    manifest,
                    _tumorNormalOptions.SomaticVcf,
                    outputVcfPath);
            return callSet;
        }
    }
}