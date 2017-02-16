using System.Collections.Generic;
using System.Text;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Isas.Framework.DataTypes;
using Isas.Framework.Logging;
using Isas.Framework.Utilities;
using Isas.Framework.WorkManagement;

namespace Canvas.Wrapper
{
    /// <summary>
    /// Run Canvas to generate CNV calls:
    /// </summary>
    public class CanvasSmallPedigreeWrapper
    {
        private readonly IWorkManager _workManager;
        private readonly ILogger _logger;
        private readonly IFileLocation _canvasExe;
        private readonly ICanvasAnnotationFileProvider _annotationFileProvider;
        private readonly ICanvasSingleSampleInputCommandLineBuilder _singleSampleInputCommandLineBuilder;
        private readonly CanvasPloidyVcfCreator _canvasPloidyVcfCreator;
        private readonly IFileLocation _mono;

        public CanvasSmallPedigreeWrapper(
            IWorkManager workManager,
            ILogger logger,
            IFileLocation canvasExe, IFileLocation mono,
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
            _mono = mono;
        }


        public CanvasSmallPedigreeOutput Run(CanvasSmallPedigreeInput input, IDirectoryLocation sampleSandbox)
        {
            var pedigreeName = input.PedigreeName;
            if (!_annotationFileProvider.IsSupported(input.GenomeMetadata))
            {
                _logger.Info($"Skipping Canvas for pedigree {pedigreeName}: unsupported reference genome '{input.GenomeMetadata.Name}'");
                return null;
            }

            if (!_annotationFileProvider.CustomDbSnpVcf(input.GenomeMetadata) && input.Vcf == null)
            {
                _logger.Info($"Skipping Canvas for sample {pedigreeName}. A dbSNP VCF file was not provided and no small variant VCF file is available");
                return null;
            }

            StringBuilder commandLine = new StringBuilder("Germline-WGS");
            commandLine.Append(_singleSampleInputCommandLineBuilder.GetSingleSampleCommandLine(pedigreeName, input.Bam, input.GenomeMetadata, sampleSandbox));

            // use normal vcf by default (performance could be similar with dbSNP vcf though)
            IFileLocation bAlleleVcf = input.Vcf.VcfFile;
            if (_annotationFileProvider.CustomDbSnpVcf(input.GenomeMetadata))
            {
                bAlleleVcf = _annotationFileProvider.GetDbSnpVcf(input.GenomeMetadata);
                commandLine.Append($" --population-b-allele-vcf {bAlleleVcf.WrapWithShellQuote()}");
            }
            else
            {
                commandLine.Append($" --sample-b-allele-vcf {bAlleleVcf.WrapWithShellQuote()}");
            }

            var ploidyVcf = _canvasPloidyVcfCreator.CreatePloidyVcf(input.PloidyInfos, input.GenomeMetadata, sampleSandbox);
            if (ploidyVcf != null)
                commandLine.Append($" --ploidy-bed {ploidyVcf.VcfFile.WrapWithShellQuote()}");
            var canvasPartitionParam = $@"--commoncnvs {_annotationFileProvider.GetCanvasAnnotationFile(input.GenomeMetadata, "commoncnvs.bed").WrapWithEscapedShellQuote()}";
            var moreCustomParameters = new Dictionary<string, string>();
            moreCustomParameters["CanvasPartition"] = canvasPartitionParam;
            commandLine.Append(_singleSampleInputCommandLineBuilder.GetCustomParameters(moreCustomParameters));
            commandLine = _singleSampleInputCommandLineBuilder.MergeCustomCanvasParameters(commandLine);

            UnitOfWork singleSampleJob = new UnitOfWork()
            {
                ExecutablePath = CrossPlatform.IsThisLinux() ? _mono.FullName : _canvasExe.FullName,
                CommandLine = CrossPlatform.IsThisLinux() ? _canvasExe + " " + commandLine : commandLine.ToString(),
                LoggingFolder = _workManager.LoggingFolder.FullName,
                LoggingStub = "Canvas_" + pedigreeName,
            };
            _workManager.DoWorkSingleThread(singleSampleJob);
            return GetCanvasOutput(pedigreeName, sampleSandbox);
        }

        private CanvasSmallPedigreeOutput GetCanvasOutput(string pedigreeName, IDirectoryLocation sampleSandbox)
        {
            var cnvVcf = new Vcf(sampleSandbox.GetFileLocation("CNV.vcf.gz"));
            var tempCnvDirectory = sampleSandbox.GetDirectoryLocation($"TempCNV_{pedigreeName}");
            var variantFrequencies = tempCnvDirectory.GetFileLocation($"VFResults{pedigreeName}.txt.gz");
            var variantFrequenciesBaf = tempCnvDirectory.GetFileLocation($"VFResults{pedigreeName}.txt.gz.baf");
            IFileLocation coverageAndVariantFrequencies = sampleSandbox.GetFileLocation("CNV.CoverageAndVariantFrequency.txt");
            IFileLocation tempStub = tempCnvDirectory.GetFileLocation($"{pedigreeName}");
            IFileLocation partitioned = tempStub.AppendName(".partitioned");
            return new CanvasSmallPedigreeOutput(cnvVcf, coverageAndVariantFrequencies, variantFrequencies,
                variantFrequenciesBaf, partitioned);
        }
    }
}