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
            List<CanvasCallset> callSets = new List<CanvasCallset>();
            for (int i = 0; i < SmallPedigreeOptions.Samples.Count(); i++)
            {
                string sampleName = SmallPedigreeOptions.Samples[i].SampleName;
                IDirectoryLocation outputDirectory = new DirectoryLocation(Path.Combine(CommonOptions.OutputDirectory.FullName, sampleName));
                Directory.CreateDirectory(outputDirectory.FullName);
                IFileLocation outputVcfPath = outputDirectory.GetFileLocation("CNV.vcf.gz");
                CanvasCallset callSet = new CanvasCallset(
                    SmallPedigreeOptions.Samples[i].Bam,
                    sampleName,
                    CommonOptions.WholeGenomeFasta,
                    CommonOptions.OutputDirectory,
                    CommonOptions.KmerFasta,
                    CommonOptions.FilterBed,
                    SmallPedigreeOptions.Samples[i].PloidyVcf,
                    SmallPedigreeOptions.Samples[i].BAlleleSites,
                    CommonOptions.IsDbSnpVcf,
                    Enumerable.Empty<IFileLocation>(),
                    null,
                    null,
                    outputVcfPath);
                callSets.Add(callSet);
            }
            return new SmallPedigreeCallset(callSets, SmallPedigreeOptions.CommonCnvsBed, SmallPedigreeOptions.PedigreeInfo);
        }
    }
}
