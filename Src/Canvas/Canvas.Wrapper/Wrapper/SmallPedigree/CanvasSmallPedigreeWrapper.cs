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
        private readonly IFileLocation _runtimeExecutable;

        public CanvasSmallPedigreeWrapper(
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

        public StringBuilder GetMultiSampleCommandLine(SampleSet<CanvasPedigreeSample> samples, GenomeMetadata genomeMetadata, Vcf vcf, IDirectoryLocation sampleSandbox)
        {
            StringBuilder commandLine = new StringBuilder();
            foreach (var sampleKvp in samples)
            {
                var sampleId = sampleKvp.Key.Id;
                var sample = sampleKvp.Value;
                commandLine.Append($" --bam \"{sample.Bam.BamFile}\" {sample.SampleType} {sampleId}");
            }
            IFileLocation kmerFasta = _annotationFileProvider.GetKmerFasta(genomeMetadata);
            commandLine.Append($" --reference \"{kmerFasta}\"");
            IDirectoryLocation wholeGenomeFasta = new FileLocation(genomeMetadata.Sequences.First().FastaPath).Directory;
            commandLine.Append($" --genome-folder \"{wholeGenomeFasta}\"");
            IFileLocation filterBed = _annotationFileProvider.GetFilterBed(genomeMetadata);
            commandLine.Append($" --filter-bed \"{filterBed}\"");
            commandLine.Append($" --output \"{sampleSandbox}\"");
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

            // use sample vcf by default (performance could be similar with dbSNP vcf though)
            var bAlleleVcf = input.Vcf.VcfFile;
            var bAlleleVcfOptionName = SingleSampleCommonOptionsParser.SampleBAlleleVcfOptionName;
            if (_annotationFileProvider.CustomDbSnpVcf(input.GenomeMetadata))
            {
                bAlleleVcf = _annotationFileProvider.GetDbSnpVcf(input.GenomeMetadata);
                bAlleleVcfOptionName = SingleSampleCommonOptionsParser.PopulationBAlleleVcfOptionName;
            }
            commandLine.Append($" --{bAlleleVcfOptionName} \"{bAlleleVcf}\"");

            var ploidyInfos = input.Samples.SelectData(sample => sample.PloidyInfo);

            var ploidyVcf = _canvasPloidyVcfCreator.CreatePloidyVcf(ploidyInfos, input.GenomeMetadata, sampleSandbox);
            if (ploidyVcf != null)
                commandLine.Append($" --{SmallPedigreeOptionsParser.PloidyVcfOptionName} \"{ploidyVcf.VcfFile}\"");
            var canvasPartitionParam = $@"--commoncnvs {_annotationFileProvider.GetCanvasAnnotationFile(input.GenomeMetadata, "commoncnvs.bed")}";

            var moreCustomParameters = new Dictionary<string, string> { ["CanvasPartition"] = canvasPartitionParam };
            commandLine.Append(_singleSampleInputCommandLineBuilder.GetCustomParameters(moreCustomParameters));
            commandLine = _singleSampleInputCommandLineBuilder.MergeCustomCanvasParameters(commandLine);
            // use Proband or, when proband is not available, first sample as pedigree id
            var pedigreeId = input.Samples.Where(x => x.Value.SampleType == SampleType.Proband).Select(x => x.Key.Id).FirstOrDefault();
            if (pedigreeId.IsNullOrEmpty())
                pedigreeId = input.Samples.First().Key.Id;

            var singleSampleJob = new UnitOfWork()
            {
                ExecutablePath = CrossPlatform.IsThisLinux() ? _runtimeExecutable.FullName : _canvasExe.FullName,
                CommandLine = CrossPlatform.IsThisLinux() ? _canvasExe + " " + commandLine : commandLine.ToString(),
                LoggingFolder = _workManager.LoggingFolder.FullName,
                LoggingStub = "Canvas_" + pedigreeId,
            };
            _workManager.DoWorkSingleThread(singleSampleJob);
            return GetCanvasOutput(input.Samples, sampleSandbox);
        }

        private CanvasSmallPedigreeOutput GetCanvasOutput(SampleSet<CanvasPedigreeSample> pedigreeSamples, IDirectoryLocation sampleSandbox)
        {
            var intermediateResults = pedigreeSamples.SelectSamples(sampleInfo =>
            {
                var sampleId = sampleInfo.Id;
                var variantFrequencies = SingleSampleCallset.GetVfSummaryPath(sampleSandbox, sampleId);
                var variantFrequenciesBaf = SingleSampleCallset.GetVfSummaryBafPath(sampleSandbox, sampleId);
                var partitioned = SingleSampleCallset.GetPartitionedPath(sampleSandbox, sampleId);
                var coverageAndVariantFrequencies = SingleSampleCallset.GetCoverageAndVariantFrequencyOutput(sampleSandbox, sampleId);
                var singleSampleVcf = SingleSampleCallset.GetSingleSamplePedigreeVcfOutput(sampleSandbox, sampleId);
                return new IntermediateOutput(new Vcf(singleSampleVcf), coverageAndVariantFrequencies, variantFrequencies, variantFrequenciesBaf, partitioned);
            });
            var cnvVcf = new Vcf(sampleSandbox.GetFileLocation("CNV.vcf.gz"));
            return new CanvasSmallPedigreeOutput(cnvVcf, intermediateResults);
        }
    }
}