using System.Linq;
using Canvas.CommandLineParsing;
using Canvas.SmallPedigree;
using CanvasCommon;
using Illumina.Common.FileSystem;

namespace Canvas
{
    public class GermlineWgsRunner : IModeRunner
    {
        private readonly GermlineWgsInput _input;

        public GermlineWgsRunner(GermlineWgsInput input)
        {
            _input = input;
        }

        public void Run(CanvasRunnerFactory runnerFactory)
        {
            var runner = runnerFactory.Create(false, CanvasCoverageMode.TruncatedDynamicRange, 100, _input.CommonOptions.CustomParams);
            var callset = GetCallset();
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