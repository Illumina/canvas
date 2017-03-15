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
            ISampleSettings sampleSettings,
            IWorkManager workManager,
            ILogger logger,
            ExecutableProcessor executableProcessor,
            DbSnpVcfProcessor dbSnpVcfProcessor,
            bool detectCnvDefault, IFileLocation pedigreefileNameStub)
        {
            _workManager = workManager;
            _logger = logger;
            _executableProcessor = executableProcessor;
            _pedigreefileNameStub = pedigreefileNameStub;
            _canvasWorkerFactory = new CanvasWorkerFactory(sampleSettings, workManager, logger, executableProcessor, dbSnpVcfProcessor, detectCnvDefault);
        }


        public ISmallPedigreeCheckpoint GetSmallPedigreeCheckpoint()
        {
            if (!_canvasWorkerFactory.RunCnvDetection()) return new NullSmallPedigreeCheckpoint(_logger);

            CanvasSmallPedigreeWrapper wrapper = new CanvasSmallPedigreeWrapper(_workManager, _logger, _canvasWorkerFactory.GetCanvasExe(), GetRuntimeExecutable(), _canvasWorkerFactory.GetAnnotationFileProvider(), _canvasWorkerFactory.GetCanvasSingleSampleInputCommandLineBuilder(_canvasWorkerFactory.GetAnnotationFileProvider()), new CanvasPloidyVcfCreator(_canvasWorkerFactory.GetPloidyCorrector()));
            return new SmallPedigreeCheckpoint(wrapper, Move, Load);
        }

        private CanvasSmallPedigreeOutput Load(CanvasSmallPedigreeInput input)
        {
            var intermediateOutputs = input.Samples.SelectData((info, sample) =>
            {
                var stub = GetStub(info.Id);
                var coverageAndVariantFrequency = SingleSampleCallset.GetCoverageAndVariantFrequencyOutput(stub);
                if (!_canvasWorkerFactory.IncludeIntermediateResults())
                    return new IntermediateOutput(coverageAndVariantFrequency, null, null, null);

                var partitioned = SingleSampleCallset.GetPartitionedPath(stub);
                var variantFrequencies = SingleSampleCallset.GetVfSummaryPath(stub);
                var variantFrequenciesBaf = SingleSampleCallset.GetVfSummaryPath(stub);
                return new IntermediateOutput(coverageAndVariantFrequency, variantFrequencies, variantFrequenciesBaf, partitioned);
            });

            return new CanvasSmallPedigreeOutput(new Vcf(GetCnvVcf()), intermediateOutputs);
        }

        private void Move(CanvasSmallPedigreeOutput source, IFileMover fileMover)
        {
            source.CnvVcf.Move(GetCnvVcf(), fileMover.Move);
            source.IntermediateOutputs.ForEach(kvp =>
            {
                MoveIntermediateOutput(kvp.Key, kvp.Value, fileMover);
            });
        }

        private IFileLocation GetCnvVcf()
        {
            return _pedigreefileNameStub.AppendName(".CNV.vcf.gz");
        }

        private IFileLocation GetStub(string sampleId)
        {
            return _pedigreefileNameStub.AppendName($"_{sampleId}.CNV");
        }

        private void MoveIntermediateOutput(SampleInfo info, IntermediateOutput output, IFileMover fileMover)
        {
            var stub = GetStub(info.Id);
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
#if DotNetCore
            return new FileLocation(_executableProcessor.GetEnvironmentExecutablePath("dotnet"));
#else
            var mono = CrossPlatform.IsThisLinux() ? _executableProcessor.GetMonoPath() : null;
            return new FileLocation(mono);
#endif
        }
    }
}