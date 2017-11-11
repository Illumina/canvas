using System;
using System.Collections.Generic;
using Canvas.CommandLineParsing;
using CanvasCommon;
using Illumina.Common.FileSystem;
using Isas.Framework.Checkpointing;
using Isas.Framework.DataTypes;
using Isas.Framework.Logging;
using Isas.Framework.WorkManagement;
using Isas.Framework.WorkManagement.CommandBuilding;

namespace Canvas.SmallPedigree
{
    public class SmallPedigreeWgsRunner : IModeRunner
    {
        public CommonOptions CommonOptions { get; }
        public SmallPedigreeOptions SmallPedigreeOptions { get; }

        public SmallPedigreeWgsRunner(SmallPedigreeInput input )
        {
            SmallPedigreeOptions = input.SmallPedigreeOptions;
            CommonOptions = input.CommonOptions;
        }

        public void Run(ILogger logger, ICheckpointRunner checkpointRunner, IWorkManager workManager, IWorkDoer workDoer, IFileLocation runtimeExecutable, Func<string, ICommandFactory> runtimeCommandPrefix)
        {
            CanvasRunner runner = new CanvasRunner(logger, workManager, workDoer, checkpointRunner, runtimeExecutable, runtimeCommandPrefix, false, CanvasCoverageMode.TruncatedDynamicRange, 100, CommonOptions.CustomParams);
            var callset = GetCallset();
            var spwWorkflow = new SmallPedigreeWorkflow(runner);
            spwWorkflow.CallPedigree(callset);
        }

        private SmallPedigreeCallset GetCallset()
        {
            var callSets = new List<PedigreeSample>();
            var outputVcf = CommonOptions.OutputDirectory.GetFileLocation("CNV.vcf.gz");
            foreach (var sample in SmallPedigreeOptions.Samples)
            {
                string sampleName = sample.SampleName;
                SingleSampleCallset callSet = new SingleSampleCallset(
                    new Bam(sample.Bam),
                    sampleName,
                    SmallPedigreeOptions.BAlleleSites,
                    SmallPedigreeOptions.IsPopulationBAlleleSites,
                    CommonOptions.OutputDirectory,
                    outputVcf);
                callSet.SampleOutputFolder.Create();
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
