using System.Collections.Generic;
using System.IO;
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
        public CommonOptions CommonOptions { get; }

        private readonly SmallPedigreeOptions _smallPedigreeOptions;

        public SmallPedigreeWgsRunner(CommonOptions commonOptions, SmallPedigreeOptions smallPedigreeOptions)
        {
            _smallPedigreeOptions = smallPedigreeOptions;
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
            List<CanvasCallset> callSets = new List<CanvasCallset>();        
            for (int i=0; i < _smallPedigreeOptions.Bams.Count(); i++) {
                    string sampleName = _smallPedigreeOptions.SampleNames.ToList()[i];
                    IDirectoryLocation outputDirectory = new DirectoryLocation(Path.Combine(CommonOptions.OutputDirectory.FullName, sampleName));
                    IFileLocation outputVcfPath = outputDirectory.GetFileLocation("CNV.vcf.gz");
                    CanvasCallset callSet = new CanvasCallset(
                    _smallPedigreeOptions.Bams.ToList()[i],
                    sampleName,
                    CommonOptions.WholeGenomeFasta,
                    CommonOptions.OutputDirectory,
                    CommonOptions.KmerFasta,
                    CommonOptions.FilterBed,
                    _smallPedigreeOptions.PloidyBed.ToList()[i],
                    _smallPedigreeOptions.BAlleleSites.ToList()[i],
                    CommonOptions.IsDbSnpVcf,
                    Enumerable.Empty<IFileLocation>(),
                    null,
                    null,
                    outputVcfPath);
                callSets.Add(callSet);
            }
            return new SmallPedigreeCallset(callSets);
        }
    }
}
