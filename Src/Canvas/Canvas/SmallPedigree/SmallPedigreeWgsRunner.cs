using System;
using System.Collections.Generic;
using System.Linq;
using Canvas.CommandLineParsing;
using CanvasCommon;
using Illumina.Common.FileSystem;
using Isas.Framework.Checkpointing;
using Isas.Framework.DataTypes;
using Isas.Framework.Logging;
using Isas.Framework.WorkManagement;

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

        public void Run(ILogger logger, ICheckpointRunner checkpointRunner, IWorkManager workManager, IFileLocation runtimeExecutable)
        {
            CanvasRunner runner = new CanvasRunner(logger, workManager, checkpointRunner, runtimeExecutable, false, CanvasCoverageMode.TruncatedDynamicRange, 100, CommonOptions.CustomParams);
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

            // Pedigree checks
            if (callSets.Count(x => x.SampleType == SampleType.Proband) > 1)
                throw new Exception("There can only be one proband in a pedigree.");

            bool haveProband = callSets.Any(x => x.SampleType == SampleType.Proband);
            bool haveMother = callSets.Any(x => x.SampleType == SampleType.Mother);
            bool haveFather = callSets.Any(x => x.SampleType == SampleType.Father);
            bool haveSibling = callSets.Any(x => x.SampleType == SampleType.Sibling);
            bool haveOther = callSets.Any(x => x.SampleType == SampleType.Other);

            bool haveTrio = haveProband && haveMother && haveFather;
            bool haveQuad = haveProband && haveMother && haveFather && haveSibling;

            if ((haveTrio || haveQuad) && haveOther)
                throw new Exception("SampleType other with trio or quad is not currently supported");
            return new SmallPedigreeCallset(callSets, analysisDetails, haveTrio);
        }
    }
}
