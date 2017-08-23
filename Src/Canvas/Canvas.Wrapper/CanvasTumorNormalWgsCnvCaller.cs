using System.Text;
using Canvas.CommandLineParsing;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Isas.Framework.DataTypes;
using Isas.Framework.Logging;
using Isas.Framework.WorkManagement;

namespace Canvas.Wrapper
{
    /// <summary>
    /// Run Canvas T/N on WGS data to generate somatic CNV calls:
    /// </summary>
    public class CanvasTumorNormalWgsCnvCaller : ICanvasCnvCaller<CanvasTumorNormalWgsInput, CanvasOutput>
    {
        private readonly IWorkManager _workManager;
        private readonly ILogger _logger;
        private readonly IFileLocation _canvasExe;
        private readonly ICanvasAnnotationFileProvider _annotationFileProvider;
        private readonly ICanvasSingleSampleInputCommandLineBuilder _singleSampleInputCommandLineBuilder;
        private readonly CanvasPloidyVcfCreator _canvasPloidyVcfCreator;
        private readonly IFileLocation _runtimeExecutable;

        public CanvasTumorNormalWgsCnvCaller(
            IWorkManager workManager,
            ILogger logger,
            IFileLocation canvasExe, IFileLocation runtimeExecutable,
            ICanvasAnnotationFileProvider annotationFileProvider,
            ICanvasSingleSampleInputCommandLineBuilder singleSampleInputCommandLineBuilder,
            CanvasPloidyVcfCreator canvasPloidyVcfCreator)
        {
            _workManager = workManager;
            _logger = logger;
            _canvasExe = canvasExe;
            _annotationFileProvider = annotationFileProvider;
            _singleSampleInputCommandLineBuilder = singleSampleInputCommandLineBuilder;
            _canvasPloidyVcfCreator = canvasPloidyVcfCreator;
            _runtimeExecutable = runtimeExecutable;
        }

        public SampleSet<CanvasOutput> Run(SampleSet<CanvasTumorNormalWgsInput> inputs, IDirectoryLocation sandbox)
        {
            var outputs = new SampleSet<CanvasOutput>();
            foreach (var input in inputs)
            {
                IDirectoryLocation sampleSandbox = sandbox.CreateSubdirectory(input.Key.Id);
                outputs.Add(input.Key, RunSingleSample(input.Key.Id, input.Value, sampleSandbox));
            }
            return outputs;
        }

        private CanvasOutput RunSingleSample(string sampleId, CanvasTumorNormalWgsInput input, IDirectoryLocation sampleSandbox)
        {
            if (!_annotationFileProvider.IsSupported(input.GenomeMetadata))
            {
                _logger.Info($"Skipping Canvas for sample {sampleId}: unsupported reference genome '{input.GenomeMetadata.Name}'");
                return null;
            }
            StringBuilder commandLine = new StringBuilder("Somatic-WGS");
            commandLine.Append(_singleSampleInputCommandLineBuilder.GetSingleSampleCommandLine(sampleId, input.TumorBam, input.GenomeMetadata, sampleSandbox));

            commandLine.Append($" --{SingleSampleCommonOptionsParser.SampleBAlleleVcfOptionName} \"{input.NormalVcf.VcfFile}\"");

            commandLine.Append($" --somatic-vcf \"{input.SomaticVcf.VcfFile}\"");

            AddPloidyVcf(commandLine, SingleSampleCommonOptionsParser.PloidyVcfOptionName, input, sampleId, sampleSandbox);

            commandLine.Append(_singleSampleInputCommandLineBuilder.GetCustomParameters());
            commandLine = _singleSampleInputCommandLineBuilder.MergeCustomCanvasParameters(commandLine);
            UnitOfWork singleSampleJob = new UnitOfWork()
            {
                ExecutablePath = CrossPlatform.IsThisLinux() ? _runtimeExecutable.FullName : _canvasExe.FullName,
                CommandLine = CrossPlatform.IsThisLinux() ? _canvasExe + " " + commandLine : commandLine.ToString(),
                LoggingStub = "Canvas_" + sampleId,
            };
            _workManager.DoWorkSingleThread(singleSampleJob);
            return GetCanvasOutput(sampleId, sampleSandbox);
        }

        private void AddPloidyVcf(StringBuilder commandLine, string optionName, CanvasTumorNormalWgsInput input, string sampleId, IDirectoryLocation sampleSandbox)
        {
            if (input.SexPloidy == null)
            {
                _logger.Warn("Sex chromosome ploidy not available. No ploidy will be provided to Canvas.");
                return;
            }
            var ploidyVcf = _canvasPloidyVcfCreator.CreatePloidyVcf(sampleId, input.SexPloidy, input.GenomeMetadata, sampleSandbox);
            commandLine.Append($" --{optionName} \"{ploidyVcf.VcfFile}\"");
        }

        private CanvasOutput GetCanvasOutput(string sampleId, IDirectoryLocation sampleSandbox)
        {
            var cnvVcf = new Vcf(sampleSandbox.GetFileLocation("CNV.vcf.gz"));
            var tempCnvDirectory = sampleSandbox.GetDirectoryLocation($"TempCNV_{sampleId}");
            var variantFrequencies = tempCnvDirectory.GetFileLocation($"VFResults{sampleId}.txt.gz");
            var variantFrequenciesBaf = tempCnvDirectory.GetFileLocation($"VFResults{sampleId}.txt.gz.baf");
            IFileLocation coverageAndVariantFrequencies = sampleSandbox.GetFileLocation("CNV.CoverageAndVariantFrequency.txt");
            IFileLocation tempStub = tempCnvDirectory.GetFileLocation($"{sampleId}");
            IFileLocation partitioned = tempStub.AppendName(".partitioned");
            return new CanvasOutput(cnvVcf, coverageAndVariantFrequencies, variantFrequencies,
                variantFrequenciesBaf, partitioned);
        }
    }
}
