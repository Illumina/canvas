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
        private readonly IWorkDoer _workDoer;
        private readonly ILogger _logger;
        private readonly IFileLocation _canvasExe;
        private readonly ICanvasAnnotationFileProvider _annotationFileProvider;
        private readonly ICanvasSingleSampleInputCommandLineBuilder _singleSampleInputCommandLineBuilder;
        private readonly CanvasPloidyVcfCreator _canvasPloidyVcfCreator;
        private readonly IFileLocation _runtimeExecutable;

        public CanvasTumorNormalWgsCnvCaller(
            IWorkDoer workDoer,
            ILogger logger,
            IFileLocation canvasExe, IFileLocation runtimeExecutable,
            ICanvasAnnotationFileProvider annotationFileProvider,
            ICanvasSingleSampleInputCommandLineBuilder singleSampleInputCommandLineBuilder,
            CanvasPloidyVcfCreator canvasPloidyVcfCreator)
        {
            _workDoer = workDoer;
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

            if (input.SexPloidy == null)
            {
                _logger.Warn("Sex chromosome ploidy not available. No ploidy will be provided to Canvas.");
            }
            else
            {

                _canvasPloidyVcfCreator.AddPloidyVcfOption(commandLine, SingleSampleCommonOptionsParser.PloidyVcfOptionName, input.GenomeMetadata, input.SexPloidy, sampleId, sampleSandbox);
            }


            commandLine.Append(_singleSampleInputCommandLineBuilder.GetCustomParameters());
            commandLine = _singleSampleInputCommandLineBuilder.MergeCustomCanvasParameters(commandLine);
            UnitOfWork singleSampleJob = new UnitOfWork()
            {
                ExecutablePath = _runtimeExecutable.FullName,
                CommandLine = _canvasExe + " " + commandLine,
                LoggingStub = "Canvas_" + sampleId,
            };
            var job = new JobInfo(_runtimeExecutable.FullName, _canvasExe + " " + commandLine, "Canvas_" + sampleId);
            return _workDoer.DoWork(WorkResourceRequest.CreateExact(1, 8), job, GetCanvasOutput(sampleId, sampleSandbox)).Await();
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
