using System.Collections.Generic;
using System.IO;
using System.Linq;
using Canvas.CommandLineParsing;
using Canvas.SmallPedigree;
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

        public SmallPedigreeOptions SmallPedigreeOptions { get; }

        public SmallPedigreeWgsRunner(CommonOptions commonOptions, SmallPedigreeOptions smallPedigreeOptions)
        {
            SmallPedigreeOptions = smallPedigreeOptions;
            CommonOptions = commonOptions;
        }

        public void Run(ILogger logger, ICheckpointRunnerAsync checkpointRunner, IWorkManager workManager)
        {
            CanvasRunner runner = new CanvasRunner(logger, workManager, checkpointRunner, false, CanvasCoverageMode.TruncatedDynamicRange, 100, CommonOptions.CustomParams);
            var callset = GetCallset();
            runner.CallPedigree(callset);
            var spwWorkflow = new SmallPedigreeWorkflow(runner);
            spwWorkflow.CallPedigree(callset);
        }

        private SmallPedigreeCallset GetCallset()
        {
            var callSets = new List<SingleSampleCallset>();
            foreach (var sample in SmallPedigreeOptions.Samples)
            {
                string sampleName = sample.SampleName;
                IDirectoryLocation outputDirectory = new DirectoryLocation(Path.Combine(CommonOptions.OutputDirectory.FullName, sampleName));
                Directory.CreateDirectory(outputDirectory.FullName);
                IFileLocation outputVcfPath = outputDirectory.GetFileLocation("CNV.vcf.gz");
                CanvasCallset callSet = new CanvasCallset(
                    sample.Bam,
                    sampleName,
                    CommonOptions.WholeGenomeFasta,
                    CommonOptions.OutputDirectory,
                    CommonOptions.KmerFasta,
                    CommonOptions.FilterBed,
                    SmallPedigreeOptions.MultiSamplePloidyVcf,
                    SmallPedigreeOptions.BAlleleSites,
                    SmallPedigreeOptions.IsPopulationBAlleleSites,
                    Enumerable.Empty<IFileLocation>(),
                    null,
                    null,
                    outputVcfPath);
                callSets.Add(new SingleSampleCallset(callSet, sample.SampleType));
            }
            return new SmallPedigreeCallset(callSets, SmallPedigreeOptions.CommonCnvsBed);
        }
    }
}
