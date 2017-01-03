using Canvas.CommandLineParsing;
using CanvasCommon;
using Illumina.SecondaryAnalysis;
using Isas.SequencingFiles;
using Isas.Shared.Checkpointing;
using Isas.Shared.Utilities;
using Isas.Shared.Utilities.FileSystem;

namespace Canvas
{
    public class TumorNormalEnrichmentRunner : IModeRunner
    {
        private readonly TumorNormalOptions _tumorNormalOptions;
        private readonly IFileLocation _normalBam;
        private readonly IFileLocation _manifest;
        public CommonOptions CommonOptions { get; }
        public SingleSampleCommonOptions SingleSampleCommonOptions { get; }

        public TumorNormalEnrichmentRunner(CommonOptions commonOptions, SingleSampleCommonOptions singleSampleCommonOptions, TumorNormalOptions tumorNormalOptions, IFileLocation normalBam, IFileLocation manifest)
        {
            _tumorNormalOptions = tumorNormalOptions;
            _normalBam = normalBam;
            _manifest = manifest;
            CommonOptions = commonOptions;
            SingleSampleCommonOptions = singleSampleCommonOptions;
        }

        public void Run(ILogger logger, ICheckpointRunnerAsync checkpointRunner, IWorkManager workManager)
        {
            CanvasRunner runner = new CanvasRunner(logger, workManager, checkpointRunner, true, CanvasCoverageMode.GCContentWeighted, 300, CommonOptions.CustomParams);
            var callset = GetCallset(logger);
            runner.CallSample(callset);
        }

        private CanvasCallset GetCallset(ILogger logger)
        {
            AnalysisDetails analysisDetails = new AnalysisDetails(CommonOptions.OutputDirectory, CommonOptions.WholeGenomeFasta, CommonOptions.KmerFasta, CommonOptions.FilterBed, SingleSampleCommonOptions.PloidyBed, null);
            IFileLocation outputVcfPath = CommonOptions.OutputDirectory.GetFileLocation("CNV.vcf.gz");
            var manifest = new NexteraManifest(_manifest.FullName, null, logger.Error);
            CanvasCallset callSet = new CanvasCallset(
                    _tumorNormalOptions.TumorBam,
                    SingleSampleCommonOptions.SampleName,
                    SingleSampleCommonOptions.BAlleleSites,
                    SingleSampleCommonOptions.IsDbSnpVcf,
                    new[] { _normalBam },
                    manifest,
                    _tumorNormalOptions.SomaticVcf,
                    outputVcfPath,
                    analysisDetails);
            return callSet;
        }
    }
}