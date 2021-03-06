using CanvasCommon;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Isas.Framework.Checkpointing;
using Isas.Framework.DataTypes;
using Isas.Framework.Logging;
using Isas.Framework.Settings;
using Isas.Framework.WorkManagement;

namespace Canvas.Wrapper.SmallPedigree
{
    public class SmallPedigreeCheckpointFactory
    {
        private readonly IWorkDoer _workDoer;
        private readonly ILogger _logger;
        private readonly ExecutableProcessor _executableProcessor;
        private readonly CanvasWorkerFactory _canvasWorkerFactory;
        private readonly IFileLocation _pedigreefileNameStub;

        public SmallPedigreeCheckpointFactory(
            IWorkDoer workDoer,
            ILogger logger,
            ExecutableProcessor executableProcessor,
            IFileLocation pedigreefileNameStub,
            CanvasWorkerFactory canvasWorkerFactory)
        {
            _workDoer = workDoer;
            _logger = logger;
            _executableProcessor = executableProcessor;
            _pedigreefileNameStub = pedigreefileNameStub;
            _canvasWorkerFactory = canvasWorkerFactory;
        }

        public ISmallPedigreeCheckpoint GetSmallPedigreeCheckpoint()
        {
            if (!_canvasWorkerFactory.RunCnvDetection()) return new NullSmallPedigreeCheckpoint(_logger);

            CanvasSmallPedigreeWrapper wrapper = new CanvasSmallPedigreeWrapper(_workDoer, _logger, _canvasWorkerFactory.GetCanvasExe(),
                GetRuntimeExecutable(), _canvasWorkerFactory.GetAnnotationFileProvider(),
                _canvasWorkerFactory.GetCanvasSingleSampleInputCommandLineBuilder(_canvasWorkerFactory.GetAnnotationFileProvider()),
                new CanvasPloidyVcfCreator(_canvasWorkerFactory.GetPloidyCorrector()));
            return new SmallPedigreeCheckpoint(wrapper, Move, Load);
        }

        private CanvasSmallPedigreeOutput Load(CanvasSmallPedigreeInput input)
        {
            var intermediateOutputs = input.Samples.SelectData((info, sample) =>
            {
                var stub = GetSingleSampleOutputStub(info);
                var coverageAndVariantFrequency = SingleSampleCallset.GetCoverageAndVariantFrequencyOutput(stub);
                var singleSampleVcf = new Vcf(SingleSampleCallset.GetVcfOutput(stub));

                var partitioned = SingleSampleCallset.GetPartitionedPath(stub);
                var variantFrequencies = SingleSampleCallset.GetVfSummaryPath(stub);
                var variantFrequenciesBaf = SingleSampleCallset.GetVfSummaryBafPath(stub);

                var coverageBigwig = SingleSampleCallset.GetCoverageBigWig(stub);
                var bAlleleBedgraph = SingleSampleCallset.GetBAlleleBedGraph(stub);
                var copyNumberBedgraph = SingleSampleCallset.GetCopyNumberBedGraph(stub);


                return new IntermediateOutput(singleSampleVcf, coverageAndVariantFrequency, variantFrequencies, variantFrequenciesBaf, partitioned, 
                    coverageBigwig, bAlleleBedgraph, copyNumberBedgraph);
            });

            return new CanvasSmallPedigreeOutput(new Vcf(GetPedigreeVcf()), intermediateOutputs);
        }

        private void Move(CanvasSmallPedigreeOutput source, IFileMover fileMover)
        {
            source.CnvVcf.Move(GetPedigreeVcf(), fileMover.Move);
            source.IntermediateOutputs.ForEach(kvp =>
            {
                MoveIntermediateOutput(kvp.Key, kvp.Value, fileMover);
            });
        }

        private IFileLocation GetPedigreeVcf()
        {
            return _pedigreefileNameStub.AppendName(".CNV.vcf.gz");
        }

        private IFileLocation GetSingleSampleOutputStub(SampleInfo sampleInfo)
        {
            return _pedigreefileNameStub.Directory.GetFileLocation($"{sampleInfo.Name}_S{sampleInfo.Number}.CNV");
        }

        private void MoveIntermediateOutput(SampleInfo info, IntermediateOutput output, IFileMover fileMover)
        {
            var stub = GetSingleSampleOutputStub(info);
            // Output:
            fileMover.Move(output.CnvVcf.VcfFile, SingleSampleCallset.GetVcfOutput(stub));
            
            // Files for visualization:
            fileMover.Move(output.CoverageBigwig, SingleSampleCallset.GetCoverageBigWig(stub));
            var targetBAlleleBedgraph = SingleSampleCallset.GetBAlleleBedGraph(stub);
            fileMover.Move(output.BAlleleBedgraph.FileLocation, targetBAlleleBedgraph.FileLocation);
            fileMover.Move(output.BAlleleBedgraph.TabixIndex, targetBAlleleBedgraph.TabixIndex);
            var targetCopyNumbedBedgraph  = SingleSampleCallset.GetCopyNumberBedGraph(stub);
            fileMover.Move(output.CopyNumberBedgraph.FileLocation, targetCopyNumbedBedgraph.FileLocation);
            fileMover.Move(output.CopyNumberBedgraph.TabixIndex, targetCopyNumbedBedgraph.TabixIndex);

            // Deprecated files:
#pragma warning disable CS0618 // Type or member is obsolete
            fileMover.Move(output.CoverageAndVariantFrequencies, SingleSampleCallset.GetCoverageAndVariantFrequencyOutput(stub)); // Used for (non-dynamic) plotting
            fileMover.Move(output.Partitioned, SingleSampleCallset.GetPartitionedPath(stub)); // used by BSVI
            fileMover.Move(output.VariantFrequencies, SingleSampleCallset.GetVfSummaryPath(stub)); // used by BSVI
            fileMover.Move(output.VariantFrequenciesBaf, SingleSampleCallset.GetVfSummaryBafPath(stub)); // used by BSVI
#pragma warning restore CS0618 // Type or member is obsolete
        }

        private IFileLocation GetRuntimeExecutable()
        {
            return new FileLocation(_executableProcessor.GetEnvironmentExecutablePath("dotnet"));
        }
    }
}