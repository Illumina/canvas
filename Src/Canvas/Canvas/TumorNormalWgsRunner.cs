using System;
using System.Linq;
using Canvas.CommandLineParsing;
using CanvasCommon;
using Illumina.Common.FileSystem;
using Isas.Framework.Checkpointing;
using Isas.Framework.Logging;
using Isas.Framework.WorkManagement;
using Isas.Framework.WorkManagement.CommandBuilding;

namespace Canvas
{
    public class TumorNormalWgsRunner : IModeRunner
    {
        private readonly TumorNormalOptions _tumorNormalOptions;
        public CommonOptions CommonOptions { get; }
        public SingleSampleCommonOptions SingleSampleCommonOptions { get; }

        public TumorNormalWgsRunner(TumorNormalWgsInput input)
        {
            _tumorNormalOptions = input.TumorNormalOptions;
            CommonOptions = input.CommonOptions;
            SingleSampleCommonOptions = input.SingleSampleCommonOptions;
        }

        public void Run(ILogger logger, ICheckpointRunner checkpointRunner, IWorkManager workManager, IWorkDoer workDoer, IFileLocation runtimeExecutable, Func<string, ICommandFactory> runtimeCommandPrefix)
        {
            CanvasRunner runner = new CanvasRunner(logger, workManager, workDoer, checkpointRunner, runtimeExecutable, runtimeCommandPrefix, true, CanvasCoverageMode.GCContentWeighted, 100, CommonOptions.CustomParams);
            var callset = GetCallset();
            runner.CallSample(callset);
        }

        private CanvasCallset GetCallset()
        {
            AnalysisDetails analysisDetails = new AnalysisDetails(CommonOptions.OutputDirectory, CommonOptions.WholeGenomeFasta, CommonOptions.KmerFasta, CommonOptions.FilterBed, SingleSampleCommonOptions.PloidyVcf, null);
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