using System;
using System.Collections.Generic;
using System.Xml.Serialization;

namespace Isas.Shared
{
	// ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
	public class SecondaryAnalysisWorkflow
	{
		#region Serializeable Types and Members
		public AmpliconWorkflowSettings AmpliconSettings = null;
		public MetagenomicsWorkflowSettings MetagenomicsSettings = null;
		public ResequencingWorkflowSettings ResequencingSettings = null;
		public WholeGenomeRnaSeqWorkflowSettings WholeGenomeRnaSeqSettings = null;
        public ZodiacRNAWorkflowSettings ZodiacRnaSettings = null;
        public RNAQuantificationWorkflowSettings RNAQuantificationSettings = null;
		public SmallRNAWorkflowSettings SmallRNASettings = null;
		public GeneralWorkflowSettings WorkflowSettings = new GeneralWorkflowSettings();
		public MtDNAWorkflowSettings MtDNASettings = null;
		public WorkflowReportSettings ReportSettings = null;
		public string Analysis = "Unknown";
		public string WorkflowVersion;
		#endregion
		// ReSharper restore InconsistentNaming

		/// <summary>
		///     Parameterless constructor - this isn't explicitly called, but it must exist in order
		///     for this type to be serializable to xml.
		/// </summary>
		public SecondaryAnalysisWorkflow()
		{
		}

		public class WorkflowReportSettings
		{
			public string AppName;
			public string AppVersion;
			public string AnalysisName;
		    public string IsUseOnly = "RUO";
		}

		// ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
		public class AmpliconWorkflowSettings
		{
            public double? AmpliconAlignerMaxSubstitutionFraction { get; set; }
			public int AmpliconAlignerMaxIndelSize { get; set; }
			public bool StitchReads { get; set; }
			public int? MaximumAllowedProbeMismatches { get; set; }
			public bool SoftClipUnalignedBases { get; set; }
			public bool CallCNV { get; set; }
			public double PseudoCount { get; set; }
			public bool GenerateAggregates { get; set; }
			public string AlignerName { get; set; }
			public string Mode { get; set; }
			// ReSharper restore InconsistentNaming
		}

		// ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
		public class GeneralWorkflowSettings
		{
			public List<string> AdapterSequences = new List<string>();
			public List<string> AdapterSequences2 = new List<string>();
			public bool TrimFastq = false; //this setting is used to trim adapters from FASTQ files if not done during BCL to FASTQ generation
			public bool UseSGE = false; // Submit executable tasks to SGE
			public bool StitchReads = false;
			public bool FindAdaptersWithIndels = false;
			public int MaximumThreadCount = 0; // 0 means no hard limit, use up to ProcessorCount threads.
			public string PumaBamModule = null;
			public string PumaVcfModule = null;
			public string PostRunCommand = null;
			public int ReadTrimQualityScoreLevel = 0; // Set >0 to enable bwa-style read trimming
			public bool ReverseComplement = false;
			public bool ClipOverlappingReadEnds = true; //used to set Isaac setting --clip-overlapping (the pairs that have read ends overlapping each other will have the lower-quality end soft-clipped)
			public VariantCallerSettings VarCallerSettings = new VariantCallerSettings();
			public SupportedVariantCallers VariantCaller = SupportedVariantCallers.GATK;
            public MetricsDeliverable Deliverable = MetricsDeliverable.Default;
            // ReSharper restore InconsistentNaming

		    public enum MetricsDeliverable
		    {
		        Default,
		        Services,
		        GEL,
		        Dx
		    }

		    public enum FastqOption
			{
				Default, // Workflow may generate fastq from bcl or use bcl directly (mostly isaac)
				UseFastq, // Like 'Default', but force fastq usage in isaac workflows
				UseExternal, // Stage from non-isas files; rewrite
				UseExisting // Stage from isas generated fastq; do not rewrite
			};
			public FastqOption FastqGeneration = FastqOption.Default; // Generate fastq from run or stage existing files
			public bool UseFastq { get { return FastqGeneration != FastqOption.Default; } }

			// Normally, we won't provide per-tile statistics during secondary analysis.
			// However, it may be useful for evaluating primary analysis results and so it can be enabled through sample sheet setting PerTileStats.
			// Sometimes this setting won't work. The problematic cases are (a) in BaseSpace, where we have a pile of FASTQs from unknown source
			// run folders, or (b) when running Isaac from FASTQ files, since Isaac throws away the 
			// original read names (including tile number) from the FASTQ file.
			public bool ComputePerTileStats = false;
			public bool CompressBam = false; // if true, provide a CRAM file in addition to a bam file. 
			public int UMILength = 0;
			[XmlIgnore]
			public Dictionary<string, string> CustomExecutablePaths = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
			[XmlIgnore]
			public Dictionary<string, string> CustomParameters = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
		}

		// ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
		public class MetagenomicsWorkflowSettings
		{
			public string TaxonomyFile = null;
			public bool DoAggregateAnalysis = false; // stores whether we should perform aggregate analysis or not
			public bool WritePerSampleStats = false; // whether we should write a stats file for each sample that duplicates info already
			// in the "MetagenomicsRunStatistics.xml" file. This is used for BaseSpace analysis.
			// ReSharper restore InconsistentNaming
		}

		// ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
		public class ResequencingWorkflowSettings
		{
            public int IsaacMajorVersion = 3;
			public SupportedAligners Aligner = SupportedAligners.bwa;
			public string ExtraBamHeader = null;
			public bool FlagPCRDuplicates = true;
			public IndelRealignment PerformIndelRealignment = IndelRealignment.None;
			public SupportedCNVCallers CNVDetection = SupportedCNVCallers.None;
			public SupportedSVCallers SVDetection = SupportedSVCallers.None;
			public BamCollapsingSettings BamCollapserSettings = new BamCollapsingSettings();
			public bool PicardHsMetrics = false;
			public string BaitManifestFileName = null;
			public int ManifestPaddingSize = 0; //default padding size is set to 150 in EnrichmentWorker
			public bool StarlingEnrichmentWithManifest = false;
			public bool FilterOffTargetVariants = true;
			public bool RunBwaAln = false;
			public bool RunExpansionHunter = false;
            public bool RunHLATyping = false;
            public ConsanguinityRunMode RunConsanguinity = ConsanguinityRunMode.ROHOnly;
			public bool GenerateAggregates = true;
            // For Neo exomes, the number of regions is over 200,000, and MSR GUI will crash if we send in such a large amount of data.
			// So, we limit the number of region stats reported to the .xml file.  The .csv files still include all regions.
			// For smaller enrichment panels like Wave and custom Neo, number of regions is far below 40000, so the max size will have no effect.
            public int EnrichmentMaxRegionStatisticsCount = 40000;
            public string GenotypeSitesVCF; // Force genotyping and no reference block compression at these sites
            // ReSharper restore InconsistentNaming

            public ResequencingWorkflowSettings()
			{
			}

			public ResequencingWorkflowSettings(SupportedAligners aligner)
			{
				Aligner = aligner;
			}

			public enum ConsanguinityRunMode
			{
				None = 0,        // don't run
				ROHOnly = 1,     // ROH calling only
				Prediction = 2,  // call ROH and predict consanguinity
				// call ROH, predict consanguinity, and output a file 
				//   of annotated SNPs for plotting.
				Extra = 3
			}
		}

		/// <summary>
		///     Tie together all data which needs to be passed to the bam collapser (used for binnining UMI reads).
		/// Currently only used by the enrichment worker.
		/// </summary>
		public class BamCollapsingSettings
		{
			public enum CombineQAlg
			{
				TakeMaxQ,
				SumQs
			};

			public int MinimumReads = 50;
			public float MaximumMismatchRate = 0.05f;
			public int MaximumBaseCallQScore = 100;
			public CombineQAlg CombineQMethod = CombineQAlg.SumQs;
		}

		/// <summary>
		///     Tie together all data which needs to be passed to the variant caller:
		/// </summary>
		// ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
		public class VariantCallerSettings
		{
			//filters only
			public VariantCallingFilterSettings Filters = new VariantCallingFilterSettings();
			public VariantCallingCombinePoolSettings HowToCombineVariantsAcrossPools = new VariantCallingCombinePoolSettings(); //Dual Strand workflows only.
			public SupportedVariantAnnotators Annotation = SupportedVariantAnnotators.None; //default set in IsasConfiguration
			public SupportedVariantAnnotators SVAnnotation = SupportedVariantAnnotators.None; //default set in IsasConfiguration
			public string TranscriptSource = "Ensembl"; // For Nirvana: "Ensembl" or "RefSeq"
			public string AnnotationDatabaseVersion; // Lock this down for a particular Isas version so results are consistent even if the service adds newer database versions. 
			public string IONAAnnotationDatabaseVersion; // Lock this down for a particular Isas version so results are consistent even if the service adds newer database versions. 
			public List<string> SkipAnnotationForVariantFilter; //annotation will be skipped for variants with these filters
			public bool UseSomaticQScoreRecalibration = false; // For Pisces only: True/false
			public float SomaticQScoreRecalibrationThreshold = 3f; // For Pisces only.  the threshold at which the recalibation alg decides to kick in.
			//how reads are handled:
			public bool StitchReads = false; // PISCES, only if XC tag present in bams
			//output only
			public bool OutputGenomeVcf = false;
			public bool GVcfBlockCompression = true;
		    public bool RetainIntermediateCNVFiles = false;
		    // ReSharper restore InconsistentNaming
		}

		// ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
		public class WholeGenomeRnaSeqWorkflowSettings
		{
			public enum Aligner { STAR, TopHat };
			public Aligner aligner = Aligner.TopHat;

			public bool FusionCalling = false;

			public bool NovelTranscripts = false;
			public bool VariantCalling = false;
			public int Trim5 = 2;
			public int Trim3 = 0;
			public int TrimToLength = 0;
			public int MinTrimmedReadLength = 24;
			public int QcMaxN = -1; // Maximum number of reads (read pairs) per sample to align
			public bool isEnrichment = false;
			public float GeneManifestOverlapCutoff = 0.7f; // Overlap cutoff used to get a gene list from intersecting a manifest with annotation
            public int? StatsSampledReads = null; // Overwrite how many reads to sample in RnaReplicateStats
			#region TopHat options
			public bool UseBowtie1 = true;
			public int? TophatReadMismatches = null;
			public int? TophatReadGapLength = null;
			public int? TophatReadEditDist = null;
			public int? TophatMateInnerDist = null;
			public int? TophatMateStdDev = null;
			public int? TophatMinIntronLength = null;
			public int? TophatMaxIntronLength = null;
			public int? TophatMaxInsertionLength = null;
			public int? TophatMaxDeletionLength = null;
			#endregion

			#region Cufflinks options
			// cuffquant/cufflinks options
			public bool? CuffquantLinksFragBiasCorrect = null;
			public bool? CuffquantLinksMultiReadCorrect = null;
			// cuffquant/cufflinks/cuffdiff options
			public bool? CuffquantLinksDiffNoEffectiveLengthCorrection = null;
			public static HashSet<string> CufflinksLibraryNormMethods = new HashSet<string>() { "classic-fpkm", "geometric", "quartile" };
			// cuffdiff/cuffnorm options
			public string CuffdiffNormLibraryNormMethod = null; // classic-fpkm, quartile
			// cuffnorm options
			public bool? CuffnormTotalHitsNorm = null;
			public bool? CuffnormCompatibleHitsNorm = null;
			// cufflinks options
			public bool? CufflinksTotalHitsNorm = null;
			public bool? CufflinksCompatibleHitsNorm = null;
			public float? CufflinksMinIsoformFraction = null;
			public float? CufflinksPreMrnaFraction = null;
			public int? CufflinksMaxIntronLength = null;
			public int? CufflinksMinIntronLength = null;
			public int? CufflinksMinFragsPerTransfrag = null;
			// cuffdiff options
			public bool? CuffdiffTotalHitsNorm = null;
			public bool? CuffdiffCompatibleHitsNorm = null;
			public float? CuffdiffFdr = null;
			public static HashSet<string> CuffdiffDispersionMethods = new HashSet<string>() { "pooled", "per-condition" };
			public string CuffdiffDispersionMethod = null; // per-condition: only available when all conditions have replicates.
			public float? CuffdiffMinIsoformFraction = null;
            #endregion

            // Additional parameterFile options passed to STAR            
            public struct StarOption { public string Option; public string Value; } // Cannot use Dict or Tuples because of serialization
            public List<StarOption> StarOptions = new List<StarOption>();

			public int ManifestPaddingSize = 0; //default padding size is set to 150

			// ReSharper restore InconsistentNaming

			protected List<DifferentialAnalysis> DifferentialAnalyses = new List<DifferentialAnalysis>();

			public List<DifferentialAnalysis> GetDifferentialAnalyses()
			{
				return DifferentialAnalyses;
			}

			public void AddDifferentialAnalysis(DifferentialAnalysis differentialAnalysis)
			{
				DifferentialAnalyses.Add(differentialAnalysis);
			}

			/// <summary>
			///     This class encapsulates a single differential analysis. A differential anlysis will have
			///     one control sample and one comparison sample
			/// </summary>
			// ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
			public class DifferentialAnalysis
			{
				public string ControlSample;
				public string ComparisonSample;
				public string Name;
				// ReSharper restore InconsistentNaming

				// Parameterless constructor for serialization
				public DifferentialAnalysis()
				{ }

				public DifferentialAnalysis(string controlSampleName, string comparisonSampleName)
				{
					Name = controlSampleName + "_vs_" + comparisonSampleName;
					ControlSample = controlSampleName;
					ComparisonSample = comparisonSampleName;
				}
			}
		}

        // ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
        public class ZodiacRNAWorkflowSettings //Todo Move this into the workflow class?
        {
            public bool STARSharedIndex = true; //Use shared mem genome index for STAR
            public int Trim5 = 0; //BP to trim from 5' and 3' end of reads before alignment
            public int Trim3 = 0;
            public int QcMaxN = -1; // Maximum number of reads (read pairs) per sample to align
            public int MinTrimmedReadLength = 24;
        }

        // ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
        public class MtDNAWorkflowSettings
		{
			public bool TinCStretchFilter = true;
			public bool MultiIndelFilter = true;
			public int MappingScore = -1;
			public string Homopolymer = "";
			public int MinimumBaseCallQuality = 30;
            public int MinimumReadCount = 10;
			public double AnalyticalThreshold = 0.25;
			public double InterpretationThreshold = 0.25;
            public int MinimumMapQuality = 12;
            public int MaxQScore = 40;
            public int MapQPenaltyWeight = 100;
            public int BaseQPenaltyWeight = 100;
            public int NoisePenaltyWeight = 400;
            public string Adapter = "";
        }

        // ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
        public class RNAQuantificationWorkflowSettings //Todo Move this into the workflow class?
        {
            public int Trim5 = 0; //BP to trim from 5' and 3' end of reads before alignment
            public int Trim3 = 0;
            public bool ExprOnly = false; // Skip steps and outputs not needed for gene quantification (counts).
            public bool Fusions = false; // Run RNA fusion detection 
            public bool SharedIndex = false; // Share STAR genome index across instances running in parallel
            public int MinTrimmedReadLength = 24;
            public List<string> ControlSamples = new List<string>(); // Differential Expression
            public List<string> ComparisonSamples = new List<string>();
            public List<string> BlockingFactors = new List<string>();
        }

        public class SmallRNAWorkflowSettings
		{
			public List<DifferentialAnalysis> DifferentialAnalyses = new List<DifferentialAnalysis>();
			public bool OutputReads = true;
			public bool NovelMiRNAs = false;
			public bool ZipNonGenomicBam = false;

			/// <summary>
			///     This class encapsulates a single differential analysis. A differential anlysis will have
			///     one control sample and one comparison sample
			/// </summary>
			// ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
			public class DifferentialAnalysis
			{
				public string Group1;
				public string Group2;
				public string Name;
				// ReSharper restore InconsistentNaming
				public DifferentialAnalysis()
				{
				}
				public DifferentialAnalysis(string controlSampleName, string comparisonSampleName)
				{
					Group1 = controlSampleName;
					Group2 = comparisonSampleName;
					Name = Group1 + "_vs_" + Group2;
				}
			}
		}
	}
}
