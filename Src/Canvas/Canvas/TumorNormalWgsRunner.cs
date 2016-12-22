using System.Linq;
using Canvas.CommandLineParsing;
using CanvasCommon;
using Illumina.SecondaryAnalysis;
using Isas.Shared.Checkpointing;
using Isas.Shared.Utilities;
using Isas.Shared.Utilities.FileSystem;

namespace Canvas
{
    public class TumorNormalWgsRunner : IModeRunner
    {
        private readonly TumorNormalOptions _tumorNormalOptions;
        public CommonOptions CommonOptions { get; }
        public SingleSampleCommonOptions SingleSampleCommonOptions { get; }

        public TumorNormalWgsRunner(CommonOptions commonOptions, SingleSampleCommonOptions singleSampleCommonOptions, TumorNormalOptions tumorNormalOptions)
        {
            _tumorNormalOptions = tumorNormalOptions;
            CommonOptions = commonOptions;
            SingleSampleCommonOptions = singleSampleCommonOptions;
        }

        public void Run(ILogger logger, ICheckpointRunnerAsync checkpointRunner, IWorkManager workManager)
        {
            CanvasRunner runner = new CanvasRunner(logger, workManager, checkpointRunner, true, CanvasCoverageMode.GCContentWeighted, 100, CommonOptions.CustomParams);
            var callset = GetCallset();
            runner.CallSample(callset);
        }

        private CanvasCallset GetCallset()
        {
            AnalysisDetails analysisDetails = new AnalysisDetails(CommonOptions.OutputDirectory, CommonOptions.WholeGenomeFasta, CommonOptions.KmerFasta, CommonOptions.FilterBed, SingleSampleCommonOptions.PloidyBed, null);
            IFileLocation outputVcfPath = CommonOptions.OutputDirectory.GetFileLocation("CNV.vcf.gz");
            CanvasCallset callSet = new CanvasCallset(
                    _tumorNormalOptions.TumorBam,
                    SingleSampleCommonOptions.SampleName,
                    SingleSampleCommonOptions.BAlleleSites,
                    SingleSampleCommonOptions.IsDbSnpVcf,
                    Enumerable.Empty<IFileLocation>(),
                    null,
                    _tumorNormalOptions.SomaticVcf,
                    outputVcfPath,
                    analysisDetails);
            return callSet;
        }
    }
}