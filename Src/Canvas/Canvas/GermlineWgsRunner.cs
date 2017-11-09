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
        private readonly GermlineWgsInput _input;

        public GermlineWgsRunner(GermlineWgsInput input)
        {
            _input = input;
        }

        public void Run(ILogger logger, ICheckpointRunner checkpointRunner, IWorkManager workManager, IFileLocation runtimeExecutable)
        {
            CanvasRunner runner = new CanvasRunner(logger, workManager, checkpointRunner, runtimeExecutable, false, CanvasCoverageMode.TruncatedDynamicRange, 100, _input.CommonOptions.CustomParams);
            var callset = GetCallset();
            logger.Info($"Normal Vcf path: {callset.SingleSampleCallset.NormalVcfPath}");
            runner.CallSample(callset);
        }

        private CanvasCallset GetCallset()
        {
            IFileLocation outputVcfPath = _input.CommonOptions.OutputDirectory.GetFileLocation("CNV.vcf.gz");
            AnalysisDetails analysisDetails = new AnalysisDetails(
                _input.CommonOptions.OutputDirectory,
                _input.CommonOptions.WholeGenomeFasta,
                _input.CommonOptions.KmerFasta,
                _input.CommonOptions.FilterBed,
                _input.SingleSampleCommonOptions.PloidyVcf,
                null);
            CanvasCallset callSet = new CanvasCallset(
                _input.Bam,
                _input.SingleSampleCommonOptions.SampleName,
                _input.SingleSampleCommonOptions.BAlleleSites,
                _input.SingleSampleCommonOptions.IsDbSnpVcf,
                Enumerable.Empty<IFileLocation>(),
                null,
                null,
                outputVcfPath,
                analysisDetails);
            return callSet;
        }
    }
}