using System.Collections.Generic;
using System.Linq;
using System.Text;
using Canvas.CommandLineParsing;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Illumina.SecondaryAnalysis.VariantCalling;
using Isas.Framework.DataTypes;
using Isas.Framework.Logging;
using Isas.Framework.Utilities;
using Isas.Framework.WorkManagement;
using Isas.SequencingFiles;

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

        public StringBuilder GetMultiSampleCommandLine(SampleSet<CanvasPedigreeSample> samples, GenomeMetadata genomeMetadata, Vcf vcf, IDirectoryLocation sampleSandbox)
        {
            StringBuilder commandLine = new StringBuilder();
            foreach (var sample in samples)
            {
                commandLine.Append($" --bam {sample.Value.Bam.BamFile.WrapWithShellQuote()}");
                commandLine.Append($" --{sample.Value.SampleType.ToString().WrapWithShellQuote()} {sample.Key.Id}");
            }
            IFileLocation kmerFasta = _annotationFileProvider.GetKmerFasta(genomeMetadata);
            commandLine.Append($" --reference {kmerFasta.WrapWithShellQuote()}");
            IDirectoryLocation wholeGenomeFasta = new FileLocation(genomeMetadata.Sequences.First().FastaPath).Directory;
            commandLine.Append($" --genome-folder {wholeGenomeFasta.WrapWithShellQuote()}");
            IFileLocation filterBed = _annotationFileProvider.GetFilterBed(genomeMetadata);
            commandLine.Append($" --filter-bed {filterBed.WrapWithShellQuote()}");
            commandLine.Append($" --output {sampleSandbox.WrapWithShellQuote()}");
            commandLine.Append(!_annotationFileProvider.CustomDbSnpVcf(genomeMetadata)
                ? $" --population-b-allele-vcf {vcf.VcfFile.WrapWithShellQuote()}"
                : $" --b-allele-vcf {vcf.VcfFile.WrapWithShellQuote()}");
            return commandLine;
        }


        public CanvasSmallPedigreeOutput Run(CanvasSmallPedigreeInput input, IDirectoryLocation sampleSandbox)
        {
            if (!_annotationFileProvider.IsSupported(input.GenomeMetadata))
            {
                _logger.Info($"Skipping Canvas: unsupported reference genome '{input.GenomeMetadata.Name}'");
                return null;
            }

            if (!_annotationFileProvider.CustomDbSnpVcf(input.GenomeMetadata) && input.Vcf == null)
            {
                _logger.Info($"Skipping Canvas. A dbSNP VCF file was not provided and no small variant VCF file is available");
                return null;
            }

            var commandLine = new StringBuilder("SmallPedigree-WGS");
            commandLine.Append(GetMultiSampleCommandLine(input.Samples, input.GenomeMetadata, input.Vcf, sampleSandbox));

            // use normal vcf by default (performance could be similar with dbSNP vcf though)
            var bAlleleVcf = input.Vcf.VcfFile;
            if (_annotationFileProvider.CustomDbSnpVcf(input.GenomeMetadata))
            {
                bAlleleVcf = _annotationFileProvider.GetDbSnpVcf(input.GenomeMetadata);
                commandLine.Append($" --population-b-allele-vcf {bAlleleVcf.WrapWithShellQuote()}");
            }
            else
            {
                commandLine.Append($" --sample-b-allele-vcf {bAlleleVcf.WrapWithShellQuote()}");
            }

            var ploidyInfos = new SampleSet<SexPloidyInfo>();
            foreach (var sample in input.Samples)
                ploidyInfos.Add(sample.Key,sample.Value.PloidyInfo);

            var ploidyVcf = _canvasPloidyVcfCreator.CreatePloidyVcf(ploidyInfos, input.GenomeMetadata, sampleSandbox);
            if (ploidyVcf != null)
                commandLine.Append($" --ploidy-bed {ploidyVcf.VcfFile.WrapWithShellQuote()}");
            var canvasPartitionParam = $@"--commoncnvs {_annotationFileProvider.GetCanvasAnnotationFile(input.GenomeMetadata, "commoncnvs.bed").WrapWithEscapedShellQuote()}";

            var moreCustomParameters = new Dictionary<string, string> {["CanvasPartition"] = canvasPartitionParam};
            commandLine.Append(_singleSampleInputCommandLineBuilder.GetCustomParameters(moreCustomParameters));
            commandLine = _singleSampleInputCommandLineBuilder.MergeCustomCanvasParameters(commandLine);
            // use Proband or, when proband is not available, first sample as pedigree id
            var pedigreeId = input.Samples.Where(x => x.Value.SampleType == SampleType.Proband).Select(x => x.Key.Id).First();
            if (pedigreeId.IsNullOrEmpty())
                pedigreeId = input.Samples.First().Key.Id;

            var singleSampleJob = new UnitOfWork()
            {
                ExecutablePath = CrossPlatform.IsThisLinux() ? _mono.FullName : _canvasExe.FullName,
                CommandLine = CrossPlatform.IsThisLinux() ? _canvasExe + " " + commandLine : commandLine.ToString(),
                LoggingFolder = _workManager.LoggingFolder.FullName,
                LoggingStub = "Canvas_" + pedigreeId,
            };
            _workManager.DoWorkSingleThread(singleSampleJob);
            var sampleNames = input.Samples.Select(x => x.Key.Id);
            return GetCanvasOutput(sampleNames, sampleSandbox);
        }

        private CanvasSmallPedigreeOutput GetCanvasOutput(IEnumerable<string> pedigreeNames, IDirectoryLocation sampleSandbox)
        {
            var coverageAndVariantFrequencies = new List<IFileLocation>();
            var variantFrequencies = new List<IFileLocation>();
            var variantFrequenciesBaf = new List<IFileLocation>();
            var partitioned = new List<IFileLocation>();
            var cnvVcf = new Vcf(sampleSandbox.GetFileLocation("CNV.vcf.gz"));
            foreach(var pedigreeName in pedigreeNames) { 
                var tempCnvDirectory = sampleSandbox.GetDirectoryLocation($"TempCNV_{pedigreeName}");
                variantFrequencies.Add(tempCnvDirectory.GetFileLocation($"VFResults{pedigreeName}.txt.gz"));
                variantFrequenciesBaf.Add(tempCnvDirectory.GetFileLocation($"VFResults{pedigreeName}.txt.gz.baf"));
                coverageAndVariantFrequencies.Add(sampleSandbox.GetFileLocation("CNV.CoverageAndVariantFrequency.txt"));
                IFileLocation tempStub = tempCnvDirectory.GetFileLocation($"{pedigreeName}");
                partitioned.Add(tempStub.AppendName(".partitioned"));
            }
            return new CanvasSmallPedigreeOutput(cnvVcf, coverageAndVariantFrequencies, variantFrequencies,
                variantFrequenciesBaf, partitioned);
        }
    }
}