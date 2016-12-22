using System.Linq;
using Canvas.CommandLineParsing;
using CanvasCommon;
using Illumina.SecondaryAnalysis;
using Isas.Shared.Checkpointing;
using Isas.Shared.Utilities;
using Isas.Shared.Utilities.FileSystem;

namespace Canvas
{
    public class GermlineWgsRunner : IModeRunner
    {
        private readonly IFileLocation _bam;
        public CommonOptions CommonOptions { get; }
        public SingleSampleCommonOptions SingleSampleCommonOptions { get; }

        public GermlineWgsRunner(CommonOptions commonOptions, SingleSampleCommonOptions singleSampleCommonOptions, IFileLocation bam)
        {
            _bam = bam;
            CommonOptions = commonOptions;
            SingleSampleCommonOptions = singleSampleCommonOptions;
        }

        public void Run(ILogger logger, ICheckpointRunnerAsync checkpointRunner, IWorkManager workManager)
        {
            CanvasRunner runner = new CanvasRunner(logger, workManager, checkpointRunner, false, CanvasCoverageMode.TruncatedDynamicRange, 100, CommonOptions.CustomParams);
            var callset = GetCallset();
            logger.Info($"Normal Vcf path: {callset.SingleSampleCallset.NormalVcfPath}");
            runner.CallSample(callset);
        }

        private CanvasCallset GetCallset()
        {
            IFileLocation outputVcfPath = CommonOptions.OutputDirectory.GetFileLocation("CNV.vcf.gz");
            AnalysisDetails analysisDetails = new AnalysisDetails(
                CommonOptions.OutputDirectory,
                CommonOptions.WholeGenomeFasta,
                CommonOptions.KmerFasta,
                CommonOptions.FilterBed,
                SingleSampleCommonOptions.PloidyBed,
                null);
            CanvasCallset callSet = new CanvasCallset(
                _bam,
                SingleSampleCommonOptions.SampleName,
                SingleSampleCommonOptions.BAlleleSites,
                SingleSampleCommonOptions.IsDbSnpVcf,
                Enumerable.Empty<IFileLocation>(),
                null,
                null,
                outputVcfPath,
                analysisDetails);
            return callSet;
        }
    }
}