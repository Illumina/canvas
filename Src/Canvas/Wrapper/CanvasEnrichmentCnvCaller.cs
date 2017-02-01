using System.Collections.Generic;
using System.Linq;
using System.Text;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Illumina.SecondaryAnalysis.VariantCalling;
using Isas.Framework.DataTypes;
using Isas.Framework.Logging;
using Isas.Framework.Utilities;
using Isas.Framework.WorkManagement;

namespace Canvas.Wrapper
{
    /// <summary>
    /// Run Canvas on enrichment data to generate CNV calls:
    /// </summary>
    public class CanvasEnrichmentCnvCaller : ICanvasCnvCaller<CanvasEnrichmentInput, CanvasEnrichmentOutput>
    {
        private readonly IWorkManager _workManager;
        private readonly ILogger _logger;
        private readonly IFileLocation _canvasExe;
        private readonly ICanvasAnnotationFileProvider _annotationFileProvider;
        private readonly ICanvasSingleSampleInputCommandLineBuilder _singleSampleInputCommandLineBuilder;
        private readonly CanvasEnrichmentInputCreator<CanvasEnrichmentInput> _enrichmentInputCreator;
        private readonly CanvasPloidyBedCreator _canvasPloidyBedCreator;

        public CanvasEnrichmentCnvCaller(IWorkManager workManager, ILogger logger, IFileLocation canvasExe,
            ICanvasAnnotationFileProvider annotationFileProvider,
            ICanvasSingleSampleInputCommandLineBuilder singleSampleInputCommandLineBuilder,
            CanvasEnrichmentInputCreator<CanvasEnrichmentInput> enrichmentInputCreator,
            CanvasPloidyBedCreator canvasPloidyBedCreator)
        {
            _workManager = workManager;
            _logger = logger;
            _canvasExe = canvasExe;
            _annotationFileProvider = annotationFileProvider;
            _singleSampleInputCommandLineBuilder = singleSampleInputCommandLineBuilder;
            _enrichmentInputCreator = enrichmentInputCreator;
            _canvasPloidyBedCreator = canvasPloidyBedCreator;
        }

        public SampleSet<CanvasEnrichmentOutput> Run(SampleSet<CanvasEnrichmentInput> inputs, IDirectoryLocation sandbox)
        {
            var targetedDbSnpVcfs = _enrichmentInputCreator.CreateDbSnpVcfForManifests(inputs, sandbox, _annotationFileProvider);
            var manifests = _enrichmentInputCreator.WriteManifests(inputs, sandbox);
            var outputs = new SampleSet<CanvasEnrichmentOutput>();
            foreach (var input in inputs)
            {
                IDirectoryLocation sampleSandbox = sandbox.CreateSubdirectory(input.Key.Id);
                var dbSnpVcf = targetedDbSnpVcfs[input.Value.NexteraManifest];
                var manifest = manifests[input.Value.NexteraManifest];
                outputs.Add(input.Key, RunSingleSample(input.Key.Id, input.Value, dbSnpVcf, manifest, sampleSandbox));
            }
            return outputs;
        }

        private CanvasEnrichmentOutput RunSingleSample(string sampleId, CanvasEnrichmentInput input, IFileLocation dbSnpVcf, IFileLocation manifest, IDirectoryLocation sampleSandbox)
        {
            if (!_annotationFileProvider.IsSupported(input.GenomeMetadata))
            {
                _logger.Info($"Skipping Canvas for sample {sampleId}: unsupported reference genome '{input.GenomeMetadata.Name}'");
                return null;
            }

            var moreCustomParameters = new Dictionary<string, string>();
            StringBuilder commandLine = new StringBuilder("Somatic-Enrichment");
            commandLine.Append(_singleSampleInputCommandLineBuilder.GetSingleSampleCommandLine(sampleId, input.Bam, input.GenomeMetadata, sampleSandbox));

            string sexChromosomeKaryotype = null;
            if (input.PloidyInfo.IsPloidyAvailable)
            {
                sexChromosomeKaryotype = input.PloidyInfo.IsXYMale.Value ? PloidyCorrector.PrettyPrintPloidy(1, 1)
                    : (input.PloidyInfo.IsXXFemale.Value ? PloidyCorrector.PrettyPrintPloidy(2, 0) : null);
            }
            string controlSexChromosomeKaryotype = null;
            if (input.IsCanvasNormalizePcaMode && input.PloidyInfo.IsPloidyAvailable)
            {
                IFileLocation pcaModelFile = input.PloidyInfo.IsXYMale.Value ? input.PcaModels.MaleModelFile
                    : (input.PloidyInfo.IsXXFemale.Value ? input.PcaModels.FemaleModelFile : null);
                if (pcaModelFile == null)
                {
                    string sampleSex = input.PloidyInfo.IsXYMale.Value ? "male"
                        : (input.PloidyInfo.IsXXFemale.Value ? "female" : "sex unknown");
                    _logger.Info($"Skipping Canvas for sample {sampleId}: PCA model file not available for {sampleSex} samples.");
                    return null;
                }
                moreCustomParameters["CanvasNormalize"] = "-m PCA";
                commandLine.Append($" --control-binned {pcaModelFile.WrapWithShellQuote()}");
                commandLine.Append($" --control-bin-size 100"); // use a dummy bin size for now
                controlSexChromosomeKaryotype = sexChromosomeKaryotype;
            }
            else if (input.PrecomputedControl?.BinnedPath != null)
            {
                commandLine.Append($" --control-binned {input.PrecomputedControl.BinnedPath.WrapWithShellQuote()}");
                commandLine.Append($" --control-bin-size {input.PrecomputedControl.BinSize}");
                controlSexChromosomeKaryotype = input.PrecomputedControl.SexChromosomeKaryotype;
            }
            else
            {
                foreach (var normalBam in input.NormalBamPaths)
                {
                    commandLine.Append($" --control-bam {normalBam.BamFile.WrapWithShellQuote()}");
                }
            }
            if (controlSexChromosomeKaryotype != null)
            {
                var controlPloidyBed = sampleSandbox.GetFileLocation("control-ploidy.bed.gz");
                if (_canvasPloidyBedCreator.GeneratePloidyBedFileFromSexChromosomeKaryotype(input.GenomeMetadata, input.GenomeMetadata.Sequences.First().FastaPath,
                    controlSexChromosomeKaryotype, controlPloidyBed.FullName, sampleSandbox.FullName))
                {
                    commandLine.Append($" --control-ploidy-bed {controlPloidyBed.WrapWithShellQuote()}");
                }
            }
            commandLine.Append($" --population-b-allele-vcf {dbSnpVcf.WrapWithShellQuote()}");
            commandLine.Append($" --manifest {manifest.WrapWithShellQuote()}");

            if (sexChromosomeKaryotype != null)
            {
                var ploidyBed = sampleSandbox.GetFileLocation("ploidy.bed.gz");
                if (_canvasPloidyBedCreator.GeneratePloidyBedFileFromSexChromosomeKaryotype(input.GenomeMetadata, input.GenomeMetadata.Sequences.First().FastaPath,
                    sexChromosomeKaryotype, ploidyBed.FullName, sampleSandbox.FullName))
                {
                    commandLine.Append($" --ploidy-bed {ploidyBed.WrapWithShellQuote()}");
                }
            }

            if (input.PredefinedBinsFile != null)
            {
                moreCustomParameters["CanvasBin"] = $"-n={input.PredefinedBinsFile.WrapWithEscapedShellQuote()}";
            }
            commandLine.Append(_singleSampleInputCommandLineBuilder.GetCustomParameters(moreCustomParameters));
            commandLine = _singleSampleInputCommandLineBuilder.MergeCustomCanvasParameters(commandLine);

            UnitOfWork singleSampleJob = new UnitOfWork()
            {
                ExecutablePath = CrossPlatform.IsThisMono() ? Utilities.GetMonoPath() : _canvasExe.FullName,
                CommandLine = CrossPlatform.IsThisMono() ? _canvasExe + " " + commandLine : commandLine.ToString(),
                LoggingFolder = _workManager.LoggingFolder.FullName,
                LoggingStub = "Canvas_" + sampleId,
            };
            _workManager.DoWorkSingleThread(singleSampleJob);
            return GetCanvasOutput(sampleId, sampleSandbox);
        }

        private CanvasEnrichmentOutput GetCanvasOutput(string sampleId, IDirectoryLocation sampleSandbox)
        {
            var cnvVcf = new Vcf(sampleSandbox.GetFileLocation("CNV.vcf.gz"));
            var tempCnvDirectory = sampleSandbox.GetDirectoryLocation($"TempCNV_{sampleId}");
            var variantFrequencies = tempCnvDirectory.GetFileLocation($"VFResults{sampleId}.txt.gz");
            var variantFrequenciesBaf = tempCnvDirectory.GetFileLocation($"VFResults{sampleId}.txt.gz.baf");
            IFileLocation coverageAndVariantFrequencies = sampleSandbox.GetFileLocation("CNV.CoverageAndVariantFrequency.txt");
            IFileLocation tempStub = tempCnvDirectory.GetFileLocation($"{sampleId}");
            IFileLocation partitioned = tempStub.AppendName(".partitioned");
            var canvasOutput = new CanvasOutput(cnvVcf, coverageAndVariantFrequencies, variantFrequencies,
                variantFrequenciesBaf, partitioned);
            IFileLocation binSize = tempStub.AppendName(".binsize");
            IFileLocation normalBinned = tempStub.AppendName(".normal.binned");
            IFileLocation unsmoothedCnd = tempStub.AppendName(".ratio.binned.cnd");
            if (!binSize.Exists) binSize = null;
            return new CanvasEnrichmentOutput(canvasOutput, binSize, normalBinned, unsmoothedCnd);
        }
    }
}