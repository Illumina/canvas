using System;
using System.Collections.Generic;
using Illumina.Common.FileSystem;
using Illumina.SecondaryAnalysis.VariantCalling;
using Isas.Framework;
using Isas.Framework.Logging;
using Isas.Framework.Settings;
using Isas.Framework.WorkManagement;
using Isas.SequencingFiles;

namespace Canvas.Wrapper
{
    public class CanvasWorkerFactory
    {
        private readonly IWorkManager _workManager;
        private readonly ISampleSettings _sampleSettings;
        private readonly ILogger _logger;
        private readonly ExecutableProcessor _executableProcessor;
        private readonly DbSnpVcfProcessor _dbSnpVcfProcessor;
        private readonly bool _detectCnvDefault;
        public static string CanvasCountPerBinSetting = "CanvasCountsPerBin";
        public static string CanvasCoverageModeSetting = "CanvasCoverageMode";
        public static string CanvasQualityScoreThresholdSetting = "CanvasQualityScoreThreshold";

        public CanvasWorkerFactory(
            ISampleSettings sampleSettings,
            IWorkManager workManager,
            ILogger logger,
            ExecutableProcessor executableProcessor,
            DbSnpVcfProcessor dbSnpVcfProcessor,
            bool detectCnvDefault)
        {
            _workManager = workManager;
            _sampleSettings = sampleSettings;
            _logger = logger;
            _executableProcessor = executableProcessor;
            _dbSnpVcfProcessor = dbSnpVcfProcessor;
            _detectCnvDefault = detectCnvDefault;
        }

        public ICanvasWorker<CanvasEnrichmentInput, CanvasEnrichmentOutput> GetCanvasEnrichmentWorker()
        {
            var mono = new FileLocation(_executableProcessor.GetMonoPath());
            var annotationProvider = GetAnnotationFileProvider();
            var canvasCnvCaller = new CanvasEnrichmentCnvCaller(
                _workManager,
                _logger,
                GetCanvasExe(),
                mono,
                annotationProvider,
                GetCanvasSingleSampleInputCommandLineBuilderWithSomaticQualityThreshold(annotationProvider),
                new CanvasEnrichmentInputCreator<CanvasEnrichmentInput>(),
                GetCanvasPloidyBedCreator());
            return GetCanvasWorker(canvasCnvCaller, CanvasEnrichmentOutput.GetFromStub);
        }

        public ICanvasWorker<CanvasResequencingInput, CanvasOutput> GetCanvasResequencingWorker(bool smallVariantCallingDisabled)
        {
            // special case: we don't do CNV calling when variant caller is disabled and we are not using a custom dbsnp vcf
            if (!CustomDbSnpVcf() && smallVariantCallingDisabled)
            {
                _logger.Info("Not running Canvas when small variant calling is disabled, unless a custom dbSNP VCF file is provided");
                if (RunCnvDetection(false))
                    throw new ArgumentException("CNV calling must be disabled when small variant calling is disabled, unless a custom dbSNP VCF file is provided");
                return new NullCanvasWorker<CanvasResequencingInput, CanvasOutput>();
            }

            var annotationProvider = GetAnnotationFileProvider();
            var mono = new FileLocation(_executableProcessor.GetMonoPath());
            var canvasCnvCaller = new CanvasResequencingCnvCaller(
                _workManager,
                _logger,
                GetCanvasExe(),
                mono,
                annotationProvider,
                GetCanvasSingleSampleInputCommandLineBuilder(annotationProvider),
                GetCanvasPloidyBedCreator());
            return GetCanvasWorker(canvasCnvCaller, CanvasOutput.GetFromStub);
        }

        private bool CustomDbSnpVcf()
        {
            return GetDbSnpVcfPath() != null;
        }

        public ICanvasWorker<CanvasTumorNormalWgsInput, CanvasOutput> GetCanvasTumorNormalWorker()
        {
            var annotationProvider = GetAnnotationFileProvider();
            var mono = new FileLocation(_executableProcessor.GetMonoPath());
            var canvasCnvCaller = new CanvasTumorNormalWgsCnvCaller(
                _workManager,
                _logger,
                GetCanvasExe(),
                mono,
                annotationProvider,
                GetCanvasSingleSampleInputCommandLineBuilderWithSomaticQualityThreshold(annotationProvider),
                GetCanvasPloidyBedCreator());
            return GetCanvasWorker(canvasCnvCaller, CanvasOutput.GetFromStub);
        }

        public ICanvasWorker<CanvasTumorNormalEnrichmentInput, CanvasOutput> GetCanvasTumorNormalEnrichmentWorker()
        {
            var annotationProvider = GetAnnotationFileProvider();
            var mono = new FileLocation(_executableProcessor.GetMonoPath());
            var canvasCnvCaller = new CanvasTumorNormalEnrichmentCnvCaller(
                _workManager,
                _logger,
                GetCanvasExe(),
                mono,
                annotationProvider,
                GetCanvasSingleSampleInputCommandLineBuilderWithSomaticQualityThreshold(annotationProvider),
                new CanvasEnrichmentInputCreator<CanvasTumorNormalEnrichmentInput>(),
                GetCanvasPloidyBedCreator());
            return GetCanvasWorker(canvasCnvCaller, CanvasOutput.GetFromStub);
        }

        private CanvasPloidyBedCreator GetCanvasPloidyBedCreator()
        {
            return new CanvasPloidyBedCreator(_logger, _workManager,
                new PloidyCorrector(_logger, _workManager, new PloidyEstimator(_logger, _workManager, _executableProcessor.GetExecutable("samtools"), false), new Isas.ClassicBioinfoTools.Tabix.TabixWrapper(_logger, _workManager, _executableProcessor.GetExecutable("tabix")), true));
        }

        public bool RequireNormalVcf()
        {
            return !CustomDbSnpVcf() && RunCnvDetection();
        }

        public bool RunCnvDetection()
        {
            return RunCnvDetection(_detectCnvDefault);
        }

        private bool RunCnvDetection(bool detectCnvDefault)
        {
            bool detectCnvs = detectCnvDefault;
            string settingString = _sampleSettings.GetStringSetting(SampleSheetUtils.RunCNVDetection, null);
            if (!string.IsNullOrEmpty(settingString))
            {
                detectCnvs = DetectCnvs(settingString);
            }
            return detectCnvs;
        }

        private static bool DetectCnvs(string name)
        {
            string tempSetting = name.ToLowerInvariant();
            switch (tempSetting)
            {
                case "0":
                case "false":
                case "none":
                    return false;
                case "1":
                case "true":
                case "canvas":
                    return true;
                default:
                    throw new Exception($"Invalid RunCNVDetection setting: {name}");
            }
        }

        private ICanvasWorker<TCanvasInput, TCanvasOutput> GetCanvasWorker<TCanvasInput, TCanvasOutput>(ICanvasCnvCaller<TCanvasInput, TCanvasOutput> canvasCnvCaller, Func<IFileLocation, bool, TCanvasOutput> getFromStub) where TCanvasInput : ICanvasCheckpointInput where TCanvasOutput : ICanvasOutput
        {
            if (!RunCnvDetection(_detectCnvDefault)) return new NullCanvasWorker<TCanvasInput, TCanvasOutput>();

            ICanvasAnnotationFileProvider annotationFileProvider = GetAnnotationFileProvider();
            bool includeIntermediateResults = _sampleSettings.GetSetting("RetainIntermediateCNVFiles", false);
            var canvasOutputNamingConventionFactory = new CanvasOutputNamingConventionFactory<TCanvasInput, TCanvasOutput>(annotationFileProvider, includeIntermediateResults, getFromStub);
            var canvasCheckpoint = new CanvasCheckpoint<TCanvasInput, TCanvasOutput>(canvasCnvCaller, canvasOutputNamingConventionFactory);
            return new CanvasWorker<TCanvasInput, TCanvasOutput>(canvasCheckpoint);
        }

        private IFileLocation GetCanvasExe()
        {
            return new FileLocation(_executableProcessor.GetExecutable("Canvas", "Canvas"));
        }

        private ICanvasAnnotationFileProvider GetAnnotationFileProvider()
        {
            return new CanvasAnnotationFileProvider(GetDbSnpVcfPath(), new ReferenceGenomeFactory());
        }

        private IFileLocation GetDbSnpVcfPath()
        {
            return _dbSnpVcfProcessor.GetDbSnpVcfPath();
        }

        private CanvasSingleSampleInputCommandLineBuilder GetCanvasSingleSampleInputCommandLineBuilder(ICanvasAnnotationFileProvider annotationFileProvider)
        {
            var allCustomParams = CommonCustomParams();
            return new CanvasSingleSampleInputCommandLineBuilder(annotationFileProvider, allCustomParams, GetCustomCanvasParameters());
        }

        private string GetCustomCanvasParameters()
        {
            return _executableProcessor.GetExecutableParameters("Canvas");
        }

        private Dictionary<string, string> CommonCustomParams()
        {
            var canvasWorkflowExecutables = new[]
            {
                "CanvasBin", "CanvasClean", "CanvasDiploidCaller", "CanvasNormalize", "CanvasPartition", "CanvasSNV", "CanvasSomaticCaller"
            };
            Dictionary<string, string> allCustomParams = new Dictionary<string, string>();
            foreach (string executable in canvasWorkflowExecutables)
            {
                string customParams = _executableProcessor.GetExecutableParameters(executable);
                if (!Illumina.Common.StringExtensions.IsNullOrWhiteSpace(customParams))
                    allCustomParams.Add(executable, customParams);
            }
            UpdateWithCanvasCountsPerBin(allCustomParams);
            UpdateWithCoverageMode(allCustomParams);
            return allCustomParams;
        }

        private CanvasSingleSampleInputCommandLineBuilder GetCanvasSingleSampleInputCommandLineBuilderWithSomaticQualityThreshold(
            ICanvasAnnotationFileProvider annotationFileProvider)
        {
            var allCustomParams = CommonCustomParams();
            UpdateWithSomaticQualityThreshold(allCustomParams);
            return new CanvasSingleSampleInputCommandLineBuilder(annotationFileProvider, allCustomParams, GetCustomCanvasParameters());
        }

        private void UpdateWithSomaticQualityThreshold(Dictionary<string, string> allCustomParams)
        {
            int? qualityScoreThreshold = _sampleSettings.GetSetting(CanvasQualityScoreThresholdSetting, (int?)null);
            if (qualityScoreThreshold.HasValue)
            {
                UpdateCustomParametersWithSetting(allCustomParams, "CanvasSomaticCaller", $" --qualitythreshold {qualityScoreThreshold.Value}");
            }
        }

        private void UpdateWithCanvasCountsPerBin(Dictionary<string, string> allCustomParams)
        {
            int? canvasCountsPerBin = _sampleSettings.GetSetting(CanvasCountPerBinSetting, (int?)null);
            if (canvasCountsPerBin.HasValue)
            {
                UpdateCustomParametersWithSetting(allCustomParams, "CanvasBin", $" -d {canvasCountsPerBin.Value}");
            }
        }

        private void UpdateWithCoverageMode(Dictionary<string, string> allCustomParams)
        {
            string canvasCoverageMode = _sampleSettings.GetSetting(CanvasCoverageModeSetting);
            if (canvasCoverageMode != null)
            {
                UpdateCustomParametersWithSetting(allCustomParams, "CanvasBin", $" -m {canvasCoverageMode}");
            }
        }

        private void UpdateCustomParametersWithSetting(Dictionary<string, string> allCustomParams, string module, string commandLineParameters)
        {
            string customParams = "";
            string tempParams;
            if (allCustomParams.TryGetValue(module, out tempParams))
            {
                customParams = tempParams;
            }
            customParams += " " + commandLineParameters;
            allCustomParams[module] = customParams;
        }
    }
}