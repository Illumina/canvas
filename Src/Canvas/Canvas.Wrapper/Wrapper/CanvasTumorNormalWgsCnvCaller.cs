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
        private readonly CanvasPloidyBedCreator _canvasPloidyBedCreator;
        private readonly IFileLocation _runtimeExecutable;

        public CanvasTumorNormalWgsCnvCaller(
            IWorkManager workManager,
            ILogger logger,
            IFileLocation canvasExe, IFileLocation runtimeExecutable,
            ICanvasAnnotationFileProvider annotationFileProvider,
            ICanvasSingleSampleInputCommandLineBuilder singleSampleInputCommandLineBuilder,
            CanvasPloidyBedCreator canvasPloidyBedCreator)
        {
            _workManager = workManager;
            _logger = logger;
            _canvasExe = canvasExe;
            _annotationFileProvider = annotationFileProvider;
            _singleSampleInputCommandLineBuilder = singleSampleInputCommandLineBuilder;
            _canvasPloidyBedCreator = canvasPloidyBedCreator;
            _runtimeExecutable = runtimeExecutable;
        }

        public SampleSet<CanvasOutput> Run(SampleSet<CanvasTumorNormalWgsInput> inputs, IDirectoryLocation sandbox)
        {
            System.Console.WriteLine("%%% Run() start");
            var outputs = new SampleSet<CanvasOutput>();
            foreach (var input in inputs)
            {
                IDirectoryLocation sampleSandbox = sandbox.CreateSubdirectory(input.Key.Id);
                outputs.Add(input.Key, RunSingleSample(input.Key.Id, input.Value, sampleSandbox));
            }
            System.Console.WriteLine("%%% Run() returning outputs");
            return outputs;
        }

        private CanvasOutput RunSingleSample(string sampleId, CanvasTumorNormalWgsInput input, IDirectoryLocation sampleSandbox)
        {
            System.Console.WriteLine("%%% RunSingleSample() start");
            if (!_annotationFileProvider.IsSupported(input.GenomeMetadata))
            {
                _logger.Info($"Skipping Canvas for sample {sampleId}: unsupported reference genome '{input.GenomeMetadata.Name}'");
                return null;
            }
            System.Console.WriteLine("%%% RunSingleSample() 2");
            StringBuilder commandLine = new StringBuilder("Somatic-WGS");
            commandLine.Append(_singleSampleInputCommandLineBuilder.GetSingleSampleCommandLine(sampleId, input.TumorBam, input.GenomeMetadata, sampleSandbox));

            commandLine.Append($" --{SingleSampleCommonOptionsParser.SampleBAlleleVcfOptionName} \"{input.NormalVcf.VcfFile}\"");

            commandLine.Append($" --somatic-vcf \"{input.SomaticVcf.VcfFile}\"");

            System.Console.WriteLine("%%% RunSingleSample() create ploidy bed");
            IFileLocation ploidyBed = _canvasPloidyBedCreator.CreatePloidyBed(input.NormalVcf, input.GenomeMetadata, sampleSandbox);
            System.Console.WriteLine("%%% RunSingleSample() ploidy bed created");
            if (ploidyBed != null)
                commandLine.Append($" --{SingleSampleCommonOptionsParser.PloidyBedOptionName} \"{ploidyBed}\"");

            commandLine.Append(_singleSampleInputCommandLineBuilder.GetCustomParameters());
            commandLine = _singleSampleInputCommandLineBuilder.MergeCustomCanvasParameters(commandLine);
            System.Console.WriteLine("%%% RunSingleSample() construct UOW");
            UnitOfWork singleSampleJob = new UnitOfWork()
            {
                ExecutablePath = CrossPlatform.IsThisLinux() ? _runtimeExecutable.FullName : _canvasExe.FullName,
                CommandLine = CrossPlatform.IsThisLinux() ? _canvasExe + " " + commandLine : commandLine.ToString(),
                LoggingStub = "Canvas_" + sampleId,
            };
            System.Console.WriteLine("%%% RunSingleSample() DWST");
            _workManager.DoWorkSingleThread(singleSampleJob);
            System.Console.WriteLine("%%% RunSingleSample() DWST complete");
            return GetCanvasOutput(sampleId, sampleSandbox);
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
