using System.Collections.Generic;
using System.Linq;
using System.Text;
using Canvas.CommandLineParsing;
using CanvasCommon;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Isas.Framework.DataTypes;
using Isas.Framework.Logging;
using Isas.Framework.WorkManagement;
using Isas.SequencingFiles;

namespace Canvas.Wrapper.SmallPedigree
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
            foreach (var sampleKvp in samples)
            {
                var sampleId = sampleKvp.Key.Id;
                var sample = sampleKvp.Value;
                commandLine.Append($" --bam {sample.Bam.BamFile.WrapWithShellQuote()}");
                if (sample.SampleType != SampleType.Other)
                    commandLine.Append($" --{sample.SampleType.GetOptionName()} {sampleId}");
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

            var ploidyInfos = input.Samples.SelectData(sample => sample.PloidyInfo);

            var ploidyVcf = _canvasPloidyVcfCreator.CreatePloidyVcf(ploidyInfos, input.GenomeMetadata, sampleSandbox);
            if (ploidyVcf != null)
                commandLine.Append($" --ploidy-bed {ploidyVcf.VcfFile.WrapWithShellQuote()}");
            var canvasPartitionParam = $@"--commoncnvs {_annotationFileProvider.GetCanvasAnnotationFile(input.GenomeMetadata, "commoncnvs.bed").WrapWithEscapedShellQuote()}";

            var moreCustomParameters = new Dictionary<string, string> { ["CanvasPartition"] = canvasPartitionParam };
            commandLine.Append(_singleSampleInputCommandLineBuilder.GetCustomParameters(moreCustomParameters));
            commandLine = _singleSampleInputCommandLineBuilder.MergeCustomCanvasParameters(commandLine);
            // use Proband or, when proband is not available, first sample as pedigree id
            var pedigreeId = input.Samples.Where(x => x.Value.SampleType == SampleType.Proband).Select(x => x.Key.Id).FirstOrDefault();
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
            var sampleBams = input.Samples.SelectData(sample => sample.Bam);
            return GetCanvasOutput(sampleBams, sampleSandbox);
        }

        private CanvasSmallPedigreeOutput GetCanvasOutput(SampleSet<Bam> pedigreeBams, IDirectoryLocation sampleSandbox)
        {
            var readGroupSamples = pedigreeBams.SelectData(GetReadGroupSample);
            var intermediateResults = readGroupSamples.SelectData(readGroupSample =>
            {
                var variantFrequencies = SingleSampleCallset.GetVfSummaryPath(sampleSandbox, readGroupSample);
                var variantFrequenciesBaf = SingleSampleCallset.GetVfSummaryBafPath(sampleSandbox, readGroupSample);
                var partitioned = SingleSampleCallset.GetPartitionedPath(sampleSandbox, readGroupSample);
                var coverageAndVariantFrequencies = SingleSampleCallset.GetCoverageAndVariantFrequencyOutput(sampleSandbox, readGroupSample);
                return new IntermediateOutput(coverageAndVariantFrequencies, variantFrequencies, variantFrequenciesBaf, partitioned);
            });
            var cnvVcf = new Vcf(sampleSandbox.GetFileLocation("CNV.vcf.gz"));
            return new CanvasSmallPedigreeOutput(cnvVcf, intermediateResults);
        }

        private static string GetReadGroupSample(Bam bam)
        {
            using (var reader = new BamReader(bam.BamFile))
            {
                return reader.GetReadGroupSample();
            }
        }
    }
}