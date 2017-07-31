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
        private readonly IWorkManager _workManager;
        private readonly ILogger _logger;
        private readonly ExecutableProcessor _executableProcessor;
        private readonly CanvasWorkerFactory _canvasWorkerFactory;
        private readonly IFileLocation _pedigreefileNameStub;

        public SmallPedigreeCheckpointFactory(
            IWorkManager workManager,
            ILogger logger,
            ExecutableProcessor executableProcessor,
            IFileLocation pedigreefileNameStub,
            CanvasWorkerFactory canvasWorkerFactory)
        {
            _workManager = workManager;
            _logger = logger;
            _executableProcessor = executableProcessor;
            _pedigreefileNameStub = pedigreefileNameStub;
            _canvasWorkerFactory = canvasWorkerFactory;
        }


        public ISmallPedigreeCheckpoint GetSmallPedigreeCheckpoint()
        {
            if (!_canvasWorkerFactory.RunCnvDetection()) return new NullSmallPedigreeCheckpoint(_logger);

            CanvasSmallPedigreeWrapper wrapper = new CanvasSmallPedigreeWrapper(_workManager, _logger, _canvasWorkerFactory.GetCanvasExe(),
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
                var singleSampleVcf = new Vcf(SingleSampleCallset.GetSingleSamplePedigreeVcfOutput(stub));
                if (!_canvasWorkerFactory.IncludeIntermediateResults())
                    return new IntermediateOutput(singleSampleVcf, coverageAndVariantFrequency, null, null, null);

                var partitioned = SingleSampleCallset.GetPartitionedPath(stub);
                var variantFrequencies = SingleSampleCallset.GetVfSummaryPath(stub);
                var variantFrequenciesBaf = SingleSampleCallset.GetVfSummaryPath(stub);
                return new IntermediateOutput(singleSampleVcf, coverageAndVariantFrequency, variantFrequencies, variantFrequenciesBaf, partitioned);
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
            fileMover.Move(output.CnvVcf.VcfFile, SingleSampleCallset.GetSingleSamplePedigreeVcfOutput(stub));
            fileMover.Move(output.CoverageAndVariantFrequencies, SingleSampleCallset.GetCoverageAndVariantFrequencyOutput(stub));
            if (_canvasWorkerFactory.IncludeIntermediateResults())
            {
                fileMover.Move(output.Partitioned, SingleSampleCallset.GetPartitionedPath(stub));
                fileMover.Move(output.VariantFrequencies, SingleSampleCallset.GetVfSummaryPath(stub));
                fileMover.Move(output.VariantFrequenciesBaf, SingleSampleCallset.GetVfSummaryPath(stub));
            }
        }

        private IFileLocation GetRuntimeExecutable()
        {
            return new FileLocation(_executableProcessor.GetEnvironmentExecutablePath("dotnet"));
        }
    }
}