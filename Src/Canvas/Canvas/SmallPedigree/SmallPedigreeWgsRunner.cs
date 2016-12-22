using System.Collections.Generic;
using System.IO;
using System.Linq;
using Canvas.CommandLineParsing;
using Canvas.SmallPedigree;
using CanvasCommon;
using Illumina.SecondaryAnalysis;
using Isas.Shared.Checkpointing;
using Isas.Shared.DataTypes;
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
            var callSets = new List<PedigreeSample>();
            foreach (var sample in SmallPedigreeOptions.Samples)
            {
                string sampleName = sample.SampleName;
                IDirectoryLocation outputDirectory = new DirectoryLocation(Path.Combine(CommonOptions.OutputDirectory.FullName, sampleName));
                Directory.CreateDirectory(outputDirectory.FullName);
                IFileLocation outputVcfPath = CommonOptions.OutputDirectory.GetFileLocation("CNV.vcf.gz");
                SingleSampleCallset callSet = new SingleSampleCallset(
                    new Bam(sample.Bam),
                    sampleName,
                    SmallPedigreeOptions.BAlleleSites,
                    SmallPedigreeOptions.IsPopulationBAlleleSites,
                    outputDirectory,
                    outputVcfPath);
                callSets.Add(new PedigreeSample(callSet, sample.SampleType));
            }
            AnalysisDetails analysisDetails = new AnalysisDetails(
                CommonOptions.OutputDirectory, 
                CommonOptions.WholeGenomeFasta, 
                CommonOptions.KmerFasta, 
                CommonOptions.FilterBed,
                SmallPedigreeOptions.MultiSamplePloidyVcf,
                SmallPedigreeOptions.CommonCnvsBed);
            return new SmallPedigreeCallset(callSets, analysisDetails);
        }
    }
}
