using System.Collections.Generic;
using System.Linq;
using Canvas.CommandLineParsing;
using CanvasCommon;
using Illumina.SecondaryAnalysis;
using Isas.Shared.Checkpointing;
using Isas.Shared.Utilities;
using Isas.Shared.Utilities.FileSystem;

namespace Canvas
{
    public class SmallPedigreeWgsRunner : IModeRunner
    {
        private readonly IEnumerable<IFileLocation> _bam;
        public CommonOptions CommonOptions { get; }

        public SmallPedigreeWgsRunner(CommonOptions commonOptions, IEnumerable<IFileLocation> bam)
        {
            _bam = bam;
            CommonOptions = commonOptions;
        }

        public void Run(ILogger logger, ICheckpointRunnerAsync checkpointRunner, IWorkManager workManager)
        {
            CanvasRunner runner = new CanvasRunner(logger, workManager, checkpointRunner, false, CanvasCoverageMode.TruncatedDynamicRange, 100, CommonOptions.CustomParams);
            var callset = GetCallset();
            runner.CallPedigree(callset);
        }

        private SmallPedigreeCallset GetCallset()
        {
            IFileLocation outputVcfPath = CommonOptions.OutputDirectory.GetFileLocation("CNV.vcf.gz");
            SmallPedigreeCallset callSet = new SmallPedigreeCallset(
                    CommonOptions.OutputDirectory,
                    _bam, 
                    new List<string>(), 
                    CommonOptions.WholeGenomeFasta,
                    CommonOptions.KmerFasta,
                    CommonOptions.FilterBed,
                    CommonOptions.PloidyBed,
                    Enumerable.Empty<IFileLocation>(),
                    outputVcfPath);
            return callSet;
        }
    }
}
