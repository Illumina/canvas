using System;
using Canvas.CommandLineParsing;
using CanvasCommon;
using Illumina.Common.FileSystem;
using Isas.Framework.Checkpointing;
using Isas.Framework.Logging;
using Isas.Framework.WorkManagement;
using Isas.Framework.WorkManagement.CommandBuilding;
using Isas.Manifests.NexteraManifest;

namespace Canvas
{
    public class TumorNormalEnrichmentRunner : IModeRunner
    {
        private readonly TumorNormalOptions _tumorNormalOptions;
        private readonly IFileLocation _normalBam;
        private readonly IFileLocation _manifest;
        public CommonOptions CommonOptions { get; }
        public SingleSampleCommonOptions SingleSampleCommonOptions { get; }

        public TumorNormalEnrichmentRunner(TumorNormalEnrichmentInput input)
        {
            _tumorNormalOptions = input.TumorNormalOptions;
            _normalBam = input.NormalBam;
            _manifest = input.Manifest;
            CommonOptions = input.CommonOptions;
            SingleSampleCommonOptions = input.SingleSampleCommonOptions;
        }

        public void Run(ILogger logger, ICheckpointRunner checkpointRunner, IWorkManager workManager, IWorkDoer workDoer, IFileLocation runtimeExecutable, Func<string, ICommandFactory> runtimeCommandPrefix)
        {
            CanvasRunner runner = new CanvasRunner(logger, workManager, workDoer, checkpointRunner, runtimeExecutable, runtimeCommandPrefix, true, CanvasCoverageMode.GCContentWeighted, 300, CommonOptions.CustomParams);
            var callset = GetCallset(logger);
            runner.CallSample(callset);
        }

        private CanvasCallset GetCallset(ILogger logger)
        {
            AnalysisDetails analysisDetails = new AnalysisDetails(CommonOptions.OutputDirectory, CommonOptions.WholeGenomeFasta, CommonOptions.KmerFasta, CommonOptions.FilterBed, SingleSampleCommonOptions.PloidyVcf, null);
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