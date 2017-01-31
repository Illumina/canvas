using System;
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
    /// Run Canvas T/N on Enrichment data to generate somatic CNV calls:
    /// </summary>
    public class CanvasTumorNormalEnrichmentCnvCaller : ICanvasCnvCaller<CanvasTumorNormalEnrichmentInput, CanvasOutput>
    {
        private readonly IWorkManager _workManager;
        private readonly ILogger _logger;
        private readonly IFileLocation _canvasExe;
        private readonly ICanvasAnnotationFileProvider _annotationFileProvider;
        private readonly ICanvasSingleSampleInputCommandLineBuilder _singleSampleInputCommandLineBuilder;
        private readonly CanvasEnrichmentInputCreator<CanvasTumorNormalEnrichmentInput> _enrichmentInputCreator;
        private readonly CanvasPloidyBedCreator _canvasPloidyBedCreator;

        public CanvasTumorNormalEnrichmentCnvCaller(
            IWorkManager workManager,
            ILogger logger,
            IFileLocation canvasExe,
            ICanvasAnnotationFileProvider annotationFileProvider,
            ICanvasSingleSampleInputCommandLineBuilder singleSampleInputCommandLineBuilder, CanvasEnrichmentInputCreator<CanvasTumorNormalEnrichmentInput> enrichmentInputCreator, CanvasPloidyBedCreator canvasPloidyBedCreator)
        {
            _workManager = workManager;
            _logger = logger;
            _canvasExe = canvasExe;
            _annotationFileProvider = annotationFileProvider;
            _singleSampleInputCommandLineBuilder = singleSampleInputCommandLineBuilder;
            _enrichmentInputCreator = enrichmentInputCreator;
            _canvasPloidyBedCreator = canvasPloidyBedCreator;
        }

        public SampleSet<CanvasOutput> Run(SampleSet<CanvasTumorNormalEnrichmentInput> inputs, IDirectoryLocation sandbox)
        {
            var bAlleleVcfs = GetBAlleleVcfs(inputs, sandbox);
            var manifests = _enrichmentInputCreator.WriteManifests(inputs, sandbox);
            var outputs = new SampleSet<CanvasOutput>();
            foreach (var input in inputs)
            {
                IDirectoryLocation sampleSandbox = sandbox.CreateSubdirectory(input.Key.Id);
                var bAlleleVcf = bAlleleVcfs[input.Key].Item1;
                var isDbSnpVcf = bAlleleVcfs[input.Key].Item2;
                var manifest = manifests[input.Value.NexteraManifest];
                outputs.Add(input.Key, RunSingleSample(input.Key.Id, input.Value, bAlleleVcf, isDbSnpVcf, manifest, sampleSandbox));
            }
            return outputs;
        }

        public SampleSet<Tuple<IFileLocation, bool>> GetBAlleleVcfs(SampleSet<CanvasTumorNormalEnrichmentInput> inputs, IDirectoryLocation sandbox)
        {
            var inputsUsingDbSnpVcf = inputs.WhereData(input => _annotationFileProvider.CustomDbSnpVcf(input.GenomeMetadata));
            var targetedDbSnpVcfs = _enrichmentInputCreator.CreateDbSnpVcfForManifests(inputsUsingDbSnpVcf, sandbox, _annotationFileProvider);
            var bAlleleVcfs = new SampleSet<Tuple<IFileLocation, bool>>();
            foreach (var input in inputs)
            {
                IFileLocation bAlleleVcf;
                bool isDbSnpVcf = true;
                if (!targetedDbSnpVcfs.TryGetValue(input.Value.NexteraManifest, out bAlleleVcf))
                {
                    bAlleleVcf = input.Value.NormalVcf.VcfFile;
                    isDbSnpVcf = false;
                }
                bAlleleVcfs[input.Key] = Tuple.Create(bAlleleVcf, isDbSnpVcf);
            }
            return bAlleleVcfs;
        }

        private CanvasOutput RunSingleSample(string sampleId, CanvasTumorNormalEnrichmentInput input, IFileLocation bAlleleVcf, bool isDbSnpVcf, IFileLocation manifest, IDirectoryLocation sampleSandbox)
        {
            if (!_annotationFileProvider.IsSupported(input.GenomeMetadata))
            {
                _logger.Info($"Skipping Canvas for sample {sampleId}: unsupported reference genome '{input.GenomeMetadata.Name}'");
                return null;
            }

            StringBuilder commandLine = new StringBuilder("Tumor-normal-enrichment");
            commandLine.Append(_singleSampleInputCommandLineBuilder.GetSingleSampleCommandLine(sampleId, input.TumorBam, input.GenomeMetadata, sampleSandbox));
            commandLine.Append($" --normal-bam {input.NormalBam.BamFile.WrapWithShellQuote()}");
            commandLine.Append($" --manifest {manifest.WrapWithShellQuote()}");
            commandLine.Append($" --b-allele-vcf {bAlleleVcf.WrapWithShellQuote()}");
            if (isDbSnpVcf)
            {
                commandLine.Append(" --exclude-non-het-b-allele-sites");
            }

            commandLine.Append($" --somatic-vcf {input.SomaticVcf.VcfFile.WrapWithShellQuote()}");

            IFileLocation ploidyBed = _canvasPloidyBedCreator.CreatePloidyBed(input.NormalVcf, input.GenomeMetadata, sampleSandbox);
            if (ploidyBed != null)
                commandLine.Append($" --ploidy-bed {ploidyBed.WrapWithShellQuote()}");
            commandLine.Append(_singleSampleInputCommandLineBuilder.GetCustomParameters());
            commandLine = _singleSampleInputCommandLineBuilder.MergeCustomCanvasParameters(commandLine);
            UnitOfWork singleSampleJob = new UnitOfWork
            {
                ExecutablePath = CrossPlatform.IsThisMono() ? Utilities.GetMonoPath() : _canvasExe.FullName,
                CommandLine = CrossPlatform.IsThisMono() ? _canvasExe + " " + commandLine : commandLine.ToString(),
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
