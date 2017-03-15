using System.Collections.Generic;
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
    /// Run Canvas to generate CNV calls:
    /// </summary>
    public class CanvasResequencingCnvCaller : ICanvasCnvCaller<CanvasResequencingInput, CanvasOutput>
    {
        private readonly IWorkManager _workManager;
        private readonly ILogger _logger;
        private readonly IFileLocation _canvasExe;
        private readonly ICanvasAnnotationFileProvider _annotationFileProvider;
        private readonly ICanvasSingleSampleInputCommandLineBuilder _singleSampleInputCommandLineBuilder;
        private readonly CanvasPloidyBedCreator _canvasPloidyBedCreator;
        private IFileLocation _runtimeExecutable;

        public CanvasResequencingCnvCaller(
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

        public SampleSet<CanvasOutput> Run(SampleSet<CanvasResequencingInput> inputs, IDirectoryLocation sandbox)
        {
            var outputs = new SampleSet<CanvasOutput>();
            foreach (var input in inputs)
            {
                IDirectoryLocation sampleSandbox = sandbox.CreateSubdirectory(input.Key.Id);
                outputs.Add(input.Key, RunSingleSample(input.Key.Id, input.Value, sampleSandbox));
            }
            return outputs;
        }

        private CanvasOutput RunSingleSample(string sampleId, CanvasResequencingInput input, IDirectoryLocation sampleSandbox)
        {
            if (!_annotationFileProvider.IsSupported(input.GenomeMetadata))
            {
                _logger.Info($"Skipping Canvas for sample {sampleId}: unsupported reference genome '{input.GenomeMetadata.Name}'");
                return null;
            }

            if (!_annotationFileProvider.CustomDbSnpVcf(input.GenomeMetadata) && input.Vcf == null)
            {
                _logger.Info($"Skipping Canvas for sample {sampleId}. A dbSNP VCF file was not provided and no small variant VCF file is available");
                return null;
            }

            StringBuilder commandLine = new StringBuilder("Germline-WGS");
            commandLine.Append(_singleSampleInputCommandLineBuilder.GetSingleSampleCommandLine(sampleId, input.Bam, input.GenomeMetadata, sampleSandbox));

            // use sample vcf by default (performance could be similar with dbSNP vcf though)
            var bAlleleVcf = input.Vcf.VcfFile;
            var bAlleleVcfOptionName = SingleSampleCommonOptionsParser.SampleBAlleleVcfOptionName;
            if (_annotationFileProvider.CustomDbSnpVcf(input.GenomeMetadata))
            {
                bAlleleVcf = _annotationFileProvider.GetDbSnpVcf(input.GenomeMetadata);
                bAlleleVcfOptionName = SingleSampleCommonOptionsParser.PopulationBAlleleVcfOptionName;
            }
            commandLine.Append($" --{bAlleleVcfOptionName} \"{bAlleleVcf}\"");


            IFileLocation ploidyBed = _canvasPloidyBedCreator.CreateGermlinePloidyBed(input.Vcf, input.GenomeMetadata, sampleSandbox);
            if (ploidyBed != null)
                commandLine.Append($" --{SingleSampleCommonOptionsParser.PloidyBedOptionName} \"{ploidyBed}\"");
            var canvasPartitionParam = $@"--commoncnvs {_annotationFileProvider.GetCanvasAnnotationFile(input.GenomeMetadata, "commoncnvs.bed").WrapWithEscapedShellQuote()}";
            var moreCustomParameters = new Dictionary<string, string>();
            moreCustomParameters["CanvasPartition"] = canvasPartitionParam;
            commandLine.Append(_singleSampleInputCommandLineBuilder.GetCustomParameters(moreCustomParameters));
            commandLine = _singleSampleInputCommandLineBuilder.MergeCustomCanvasParameters(commandLine);

            UnitOfWork singleSampleJob = new UnitOfWork()
            {
                ExecutablePath = CrossPlatform.IsThisLinux() ? _runtimeExecutable.FullName : _canvasExe.FullName,
                CommandLine = CrossPlatform.IsThisLinux() ? _canvasExe + " " + commandLine : commandLine.ToString(),
                LoggingFolder = _workManager.LoggingFolder.FullName,
                LoggingStub = "Canvas_" + sampleId,
            };
            _workManager.DoWorkSingleThread(singleSampleJob);
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
