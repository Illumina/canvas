using System.Linq;
using Canvas.CommandLineParsing;
using CanvasCommon;
using Illumina.Common.FileSystem;
using Isas.Framework.Checkpointing;
using Isas.Framework.Logging;
using Isas.Framework.WorkManagement;

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

        public void Run(ILogger logger, ICheckpointRunner checkpointRunner, IWorkManager workManager, IFileLocation runtimeExecutable)
        {
            CanvasRunner runner = new CanvasRunner(logger, workManager, checkpointRunner, runtimeExecutable, false, CanvasCoverageMode.TruncatedDynamicRange, 100, CommonOptions.CustomParams);
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