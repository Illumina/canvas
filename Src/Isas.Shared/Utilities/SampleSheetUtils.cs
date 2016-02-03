// ReSharper disable once CheckNamespace
namespace Isas.Shared
{
    /// <summary>
    /// This class stores sample sheet setting names in a consistent place.  Going forward, we want to use this class
    /// only to store generic settings that apply to all (or at least most) workflows.  Any workflow-specific settings 
    /// should live in the same project as the worker.
    /// </summary>
    public static class SampleSheetUtils
    {
        // ReSharper disable InconsistentNaming
        public const string Adapter = "adapter";
        public const string AdapterRead2 = "adapterread2";
        public const string Aligner = "aligner";
        public const string AnalysisName = "analysisname"; //currently used for BaseSpace reports
        public const string AnnotationDatabaseVersion = "annotationdatabaseversion";
        public const string AppName = "appname";
        public const string AppVersion = "appversion";
        public const string DatasetTypeId = "datasettypeid";
        public const string BamSummaryMethod = "bamsummarymethod";
        public const string CanvasBinSizeManifest = "canvasbinsizemanifest"; // Special: We see this setting with letters appended
        public const string CanvasControlBinnedManifest = "canvascontrolbinnedmanifest"; // Special: We see this setting with letters appended
        public const string CanvasSexChromosomeKaryotypeManifest = "canvassexchromosomekaryotypemanifest"; // Special: We see this setting with letters appended
        public const string UseSGE = "usesge";
        public const string CreateFastqForIndexReads = "createfastqforindexreads";
        public const string CompressBam = "compressbam";
        public const string ComputePerTileStats = "pertilestats";
        public const string DefaultMinTrimmedReadLength = "defaultmintrimmedreadlength"; // for adapter trimming
        public const string DisablePlotGeneration = "disableplotgeneration";
        public const string ExcludeRegionsManifest = "excluderegionsmanifest"; // Special: We see this setting with letters appended
        public const string FASTQInput = "fastqinput";
        public const string FilterOffTargetVariants = "filterofftargetvariants";
        public const string FindAdaptersWithIndels = "findadapterswithindels";
        public const string FilterOutSingleStrandVariants = "filteroutsinglestrandvariants";
        public const string FilterPCRDuplicates = "filterpcrduplicates"; // Deprecated synonym of FlagPCRDuplicates; supported for back-compatibility with old sheets
        public const string FlagPCRDuplicates = "flagpcrduplicates";
        public const string GVcfBlockCompression = "gvcfblockcompression";
        public const string IndelRealignment = "indelrealignment";
        public const string IndelRepeatFilterCutoff = "indelrepeatfiltercutoff";
        public const string IsaacSeedLength = "isaacseedlength";
        public const string ManifestName = "manifestname";
        public const string MinimumCoverageDepth = "minimumcoveragedepth";
        public const string MinQScore = "minqscore";
        public const string MaximumThreadCount = "maximumthreadcount";
        public const string MaximumMemoryGB = "maximummemorygb";
        public const string OutputGenomeVCF = "outputgenomevcf";
        public const string PostRunCommand = "postruncommand";
        public const string PumaArrayVcf = "arrayvcf";
        public const string PumaBamModule = "pumabammodule";
        public const string PumaConfig = "pumaconfig";
        public const string PumaVcfModule = "pumavcfmodule";
        public const string QualityScoreTrim = "qualityscoretrim";
        public const string Read1EndWithCycle = "read1endwithcycle";
        public const string Read2EndWithCycle = "read2endwithcycle";
        public const string Read1StartFromCycle = "read1startfromcycle";
        public const string Read2StartFromCycle = "read2startfromcycle";
        public const string Read1UMILength = "read1umilength";
        public const string Read2UMILength = "read2umilength";
        public const string RetainIntermediateCNVFiles = "retainintermediatecnvfiles";
        public const string RetainTempFiles = "retaintempfiles";
        public const string ReverseComplement = "reversecomplement";
        public const string RunBwaAln = "runbwaaln";
        public const string RunCNVDetection = "runcnvdetection";
        public const string RunSVDetection = "runsvdetection";
        public const string StarlingEnrichmentWithManifest = "starlingenrichmentwithmanifest";
        public const string StrandBiasFilter = "strandbiasfilter";
        public const string SVAnnotation = "svannotation";
        public const string TranscriptSource = "transcriptsource";
        public const string TrimFastq = "trimfastq";
        public const string VariantAnnotation = "variantannotation";
        public const string VariantCaller = "variantcaller";
        public const string VariantCallerQualityModel = "variantcallerqualitymodel";
        public const string VariantCallerOnlyUseProperPairs = "variantcalleronlyuseproperpairs";
        public const string VariantCallerReportNoCalls = "variantcallerreportnocalls";
        public const string VariantCallerQScoreRecalibrationThreshold = "variantcallersomaticqscorerecalibrationthreshold";
        public const string VariantCallerUseSomaticQScoreRecalibration = "variantcallerusesomaticqscorerecalibration";
        public const string VariantCallerStitchReads = "variantcallerstitchreads";
        public const string VariantCallMaxLengthPhasedSNP = "variantcallphasedsnp_maxlength";
        public const string VariantCallMaxGapPhasedSNP = "variantcallphasedsnp_maxgap";
        public const string VariantCallPhasedSNP = "variantcallphasedsnp";
        public const string VariantFilterQualityCutoff = "variantfilterqualitycutoff"; // For back-compatibility
        public const string VariantFrequencyEmitCutoff = "variantfrequencyemitcutoff";
        public const string VariantFrequencyFilterCutoff = "variantfrequencyfiltercutoff";
        public const string VariantMinimumGQCutoff = "variantminimumgqcutoff";
        public const string VariantMinimumGQXCutoff = "variantminimumgqxcutoff";
        public const string VariantMinimumMQCutoff = "variantminimummqcutoff";
        public const string VariantMinimumQDCutoff = "variantminimumqdcutoff";
        public const string VariantMinimumQualCutoff = "variantminimumqualcutoff";
        public const string VariantQualityEmitCutoff = "variantqualityemitcutoff";
        // ReSharper restore InconsistentNaming

        public static bool IsValidSampleSheetSetting(string key)
        {
            string keyLower = key.ToLowerInvariant();

            if (keyLower.StartsWith(CanvasBinSizeManifest)) return true; // Special case: CanvasBinSizeManifest*
            if (keyLower.StartsWith(CanvasControlBinnedManifest)) return true; // Special case: CanvasControlBinnedManifest*
            if (keyLower.StartsWith(ExcludeRegionsManifest)) return true; // Special case: ExcludeRegionsManifest*
            if (keyLower.StartsWith("excludetiles")) return true; // Special case: excludetiles*
            if (keyLower.StartsWith("customexecutablepath")) return true; // Special case: CustomExecutablePath*
            if (keyLower.StartsWith("customparameters")) return true; // Special case: CustomParameters*            
            if (keyLower.StartsWith("wgrna_star_")) return true; // Special case, RNA-STAR custom options: wgrna_rna_* // todo avoid string duplication

            switch (keyLower)
            {
                case Adapter:
                case AdapterRead2:
                case Aligner:
                case AnalysisName:
                case AnnotationDatabaseVersion:
                case AppName:
                case AppVersion:
                case DatasetTypeId:
                case BamSummaryMethod:
                case CreateFastqForIndexReads:
                case UseSGE:
                case CompressBam:
                case ComputePerTileStats:
                case DefaultMinTrimmedReadLength:
                case DisablePlotGeneration:
                case FASTQInput:
                case FilterOffTargetVariants:
                case FindAdaptersWithIndels:
                case FilterOutSingleStrandVariants:
                case FilterPCRDuplicates:
                case FlagPCRDuplicates:
                case GVcfBlockCompression:
                case IndelRealignment:
                case IndelRepeatFilterCutoff:
                case IsaacSeedLength:
                case MinimumCoverageDepth:
                case MinQScore:
                case ManifestName:
                case MaximumThreadCount:
                case MaximumMemoryGB:
                case OutputGenomeVCF:
                case PostRunCommand:
                case PumaArrayVcf:
                case PumaBamModule:
                case PumaConfig:
                case PumaVcfModule:
                case QualityScoreTrim:
                case Read1EndWithCycle:
                case Read2EndWithCycle:
                case Read1StartFromCycle:
                case Read2StartFromCycle:
                case Read1UMILength:
                case Read2UMILength:
                case RetainIntermediateCNVFiles:
                case RetainTempFiles:
                case ReverseComplement:
                case RunBwaAln:
                case RunCNVDetection:
                case RunSVDetection:
                case StarlingEnrichmentWithManifest:
                case StrandBiasFilter:
                case SVAnnotation:
                case TranscriptSource:
                case TrimFastq:
                case VariantAnnotation:
                case VariantCaller:
                case VariantFilterQualityCutoff: // For back-compatibility
                case VariantFrequencyEmitCutoff:
                case VariantFrequencyFilterCutoff:
                case VariantMinimumGQCutoff:
                case VariantMinimumGQXCutoff:
                case VariantMinimumMQCutoff:
                case VariantMinimumQDCutoff:
                case VariantMinimumQualCutoff:
                case VariantQualityEmitCutoff:
                case VariantCallerQualityModel:
                case VariantCallPhasedSNP:
                case VariantCallMaxLengthPhasedSNP:
                case VariantCallerOnlyUseProperPairs:
                case "customindexprimermix": // Ignore this MCS option which sneaked its way into [Settings] :(
                case "customread1primermix": // Ignore this MCS option which sneaked its way into [Settings] :(
                case "customread2primermix": // Ignore this MCS option which sneaked its way into [Settings] :(
                case "percenttilestoscan": // Ignore this MCS option which sneaked its way into [Settings] :(
                case "specifictilestoscan": // Ignore this MCS option which sneaked its way into [Settings] :(
                    return true;
                default:
                    return false;
            }
        }
    }
}
