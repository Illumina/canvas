using System.Linq;
using Canvas.CommandLineParsing;
using CanvasCommon;
using Illumina.SecondaryAnalysis;
using Illumina.SecondaryAnalysis.Workflow;
using Isas.Shared;

namespace Canvas
{
    public class TumorNormalWgsRunner : IModeRunner
    {
        private readonly TumorNormalOptions _tumorNormalOptions;
        public CommonOptions CommonOptions { get; }

        public TumorNormalWgsRunner(CommonOptions commonOptions, TumorNormalOptions tumorNormalOptions)
        {
            _tumorNormalOptions = tumorNormalOptions;
            CommonOptions = commonOptions;
        }

        public void Run(ILogger logger, ICheckpointRunner checkpointRunner, IWorkManager workManager)
        {
            CanvasRunner runner = new CanvasRunner(logger, workManager, checkpointRunner, true, CanvasCoverageMode.GCContentWeighted, 100, CommonOptions.CustomParams);
            var callset = GetCallset();
            runner.CallSample(callset);
        }

        private CanvasCallset GetCallset()
        {
            IFileLocation outputVcfPath = CommonOptions.OutputDirectory.GetFileLocation("CNV.vcf.gz");
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
                    Enumerable.Empty<IFileLocation>(),
                    null,
                    _tumorNormalOptions.SomaticVcf,
                    outputVcfPath);
            return callSet;
        }
    }
}