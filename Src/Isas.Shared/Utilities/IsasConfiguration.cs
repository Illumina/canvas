using System;
using System.Collections.Generic;
using System.Configuration;
using System.Diagnostics;
using System.IO;
using System.Reflection;

namespace Isas.Shared
{
	/// <summary>
	///     IsasConfiguration stores the configuration for Isas.  Various default settings are provided, and can
	///     be overridden using App.config
	/// </summary>
	public class IsasConfiguration
	{
		#region Members (alphabetized!)
		public float AdapterTrimmingStringency { get; set; }
		public string BamSummaryMethod { get; set; }// "BamStats" by default. Can be set to "PUMA". If set to "Both" then BamStats will run first and PUMA will overwrite/supplement with its stats 
		public int BamSummaryThreads { get; set; }// Set to 0 to use default (6 for stats by threads with high mem or all processors for stats by processors with high mem. low mem will be 2)
		public int ClusterCoresPerSlot { get; set; }
		public float ClusterGigabytesPerSlot { get; set; }
		public string ClusterParallelEnvironmentName { get; set; }
		public string ClusterQueueName { get; set; }
		// ReSharper disable InconsistentNaming
		public bool ConvertMissingBCLsToNoCalls { get; set; }
		// ReSharper restore InconsistentNaming
		public bool CreateFastqForIndexReads { get; set; }
		public SupportedAligners DefaultAligner { get; set; }
		public SupportedVariantCallers DefaultVariantCaller { get; set; }
		public bool DisablePlotGeneration { get; set; } //can be overridden by the samplesheet setting DisablePlotGeneration
		public bool EnableBamStatsByProcess { get; set; } // true by default, can be overridden in some workflows.  one process is launched per thread. total number of threads can be specified with BamSummaryThreads 
		public bool FilterNonPFReads { get { return true; } set { if (!value) throw new ArgumentException("FilterNonPFReads must be true"); } } // Per discussion with Steven Tanner (2015-7-9) we will no longer support keeping non-PF reads.
		public int GATKDownsampleDepth { get; set; }
		public string GenomePath { get; set; }
		public int IndelRepeatFilterCutoff { get; set; }
        public int IsaacMemoryLimitFudgeFactor { get; set; }
        public float MaximumMemoryGB { get; set; }
		public float MaximumHoursPerProcess { get; set; }
		public long MaximumMegabasesIsaacIndex { get; set; }
		public int MinimumAlignReadLength { get; set; }
		public int NMaskShortAdapterReads { get; set; }
		public string PostRunCommand { get; set; }
		public bool ReportDemuxToInterops { get; set; }
		public bool RequirePerfectMiRNAHit { get; set; }
		public bool RetainTempFiles { get; set; }
		public bool ReportIndelRepeatCount { get; set; }
		public bool RunGATKForSingleSamples { get; set; }
		public SupportedVariantAnnotators SVAnnotation { get; set; }
		public string TempFolder { get; set; } // If the genome folder isn't writable, we write genomes here.
		public bool UseCluster { get; set; }
		public bool UsePipedBwa { get; set; }
		public SupportedVariantAnnotators VariantAnnotation { get; set; }
		public int VariantFilterQualityCutoff { get; set; }
		#endregion

		private static readonly object ConfigLock = new object();
		private static IsasConfiguration Config;
		private Dictionary<string, string> _dictionaryConfig;

		/// <summary>
		///     This constructor is private, to enforce a singleton IsasConfiguration object.
		///     Access it by calling GetConfiguration.
		/// </summary>
		private IsasConfiguration(Dictionary<string, string> settings)
		{
			_dictionaryConfig = settings;

			#region Set defaults (alphabetical order):
			AdapterTrimmingStringency = 0.9f; // In 2.0, this was effectively 2/3
			BamSummaryMethod = "BamStats"; //the default until all workflows have switched to PUMA. Workflow can override this by invoking UsePumaAsDefaultBamSummaryMethod.
			BamSummaryThreads = 0;
			ClusterCoresPerSlot = 1;
			ClusterGigabytesPerSlot = 4.0f;
			ClusterParallelEnvironmentName = null;
			ClusterQueueName = null;
			ConvertMissingBCLsToNoCalls = true; // *enabled* by default
			DefaultAligner = SupportedAligners.bwa;
			DefaultVariantCaller = SupportedVariantCallers.GATK;
			DisablePlotGeneration = Utilities.IsThisMono();
			EnableBamStatsByProcess = true;
			FilterNonPFReads = true;
			GATKDownsampleDepth = 5000;
            IsaacMemoryLimitFudgeFactor = 1; 
            IndelRepeatFilterCutoff = -1; // Later override this -1 with default of 8 for MiSeq, 0 for HiSeq; sample sheet overrides the instrument default
			MaximumMemoryGB = 6.5f; // Overwritten based on machine physical memory below (and later maybe overwritten by SampleSheet) 
			MaximumHoursPerProcess = 2.5f;
			MaximumMegabasesIsaacIndex = 25;
			MinimumAlignReadLength = 21;
			NMaskShortAdapterReads = 10;
			ReportDemuxToInterops = true;
			RequirePerfectMiRNAHit = true;
			RetainTempFiles = false;
			ReportIndelRepeatCount = false;
			RunGATKForSingleSamples = true;
			SVAnnotation = SupportedVariantAnnotators.Nirvana; //fall back to annotation with genes file for now
			UseCluster = false;
			UsePipedBwa = false; // Update 2/5/13: Disable piped bwa under both Windows and Linux by default, for maximum robustness.  
			VariantAnnotation = SupportedVariantAnnotators.Nirvana; // Default - Nirvana.  (No more IONA, IAS or MARS)

			// Default temp folder is "Temp" subfolder of the application folder:
			TempFolder = Path.Combine(Utilities.GetAssemblyFolder(typeof(IsasConfiguration)), "Temp");
			VariantFilterQualityCutoff = -1; //override this after we know which varcaller we are using

			#endregion

			// Attempt to supply a good default for MaximumGigabytesPerProcess using the amount of physical memory:
			// Todo: This could look at SGE environment variables to figure out the amount of RAM Isas it allowed to use
			try
			{
				double physicalGB = MachineInfo.TotalPhysicalMemoryGB();
				MaximumMemoryGB = (float)physicalGB - 2;
				// Special case: On very low-memory machines, aim a little closer to the ceiling:
				if (physicalGB < 8)
				{
					MaximumMemoryGB = (float)physicalGB - 1;
				}
			}
			catch
			{
			} // Swallow the exception - nothing to be done but fall back to the config file

			if (Utilities.IsThisMono())
			{
				//VariantAnnotation = SupportedVariantAnnotators.None; // 2/8/13: Enable annotation by default
				DefaultAligner = SupportedAligners.Isaac;
				DefaultVariantCaller = SupportedVariantCallers.Starling;
			}

			// Parse config file values!  Parse (most of them) in alphabetical order, for simplicity of maintenance:
			if (_dictionaryConfig.ContainsKey("AdapterTrimmingStringency"))
			{
				string setting = ParseSettingValue("AdapterTrimmingStringency");
				AdapterTrimmingStringency = float.Parse(setting);
				if (AdapterTrimmingStringency > 1 || AdapterTrimmingStringency < 0.5f)
				{
					throw new Exception(
						"Error: Invalid setting for AdapterTrimmingStringency - provide a value between 0.5 and 1");
				}
			}

			if (_dictionaryConfig.ContainsKey("BamSummaryMethod"))
			{
				BamSummaryMethod = ParseSettingValue("BamSummaryMethod");
			}
			if (_dictionaryConfig.ContainsKey("BamSummaryThreads"))
			{
				string setting = ParseSettingValue("BamSummaryThreads");
				BamSummaryThreads = int.Parse(setting);
			}

			if (_dictionaryConfig.ContainsKey("ClusterParallelEnvironmentName"))
			{
				ClusterParallelEnvironmentName = ParseSettingValue("ClusterParallelEnvironmentName");
			}
			if (_dictionaryConfig.ContainsKey("ClusterQueueName"))
			{
				ClusterQueueName = ParseSettingValue("ClusterQueueName");
			}
			if (_dictionaryConfig.ContainsKey("ClusterCoresPerSlot"))
			{
				ClusterCoresPerSlot = int.Parse(ParseSettingValue("ClusterCoresPerSlot"));
			}
			if (_dictionaryConfig.ContainsKey("ClusterGigabytesPerSlot"))
			{
				ClusterGigabytesPerSlot = float.Parse(ParseSettingValue("ClusterGigabytesPerSlot"));
			}

			ConvertMissingBCLsToNoCalls = GetBooleanSetting("ConvertMissingBCLsToNoCalls", ConvertMissingBCLsToNoCalls);

			if (_dictionaryConfig.ContainsKey("DefaultAligner"))
			{
				try
				{
					DefaultAligner = (SupportedAligners)
						Enum.Parse(typeof(SupportedAligners), ParseSettingValue("DefaultAligner"), true);
				}
				catch (Exception)
				{
					// not a valid setting, ignore
				}
			}

			if (_dictionaryConfig.ContainsKey("DefaultVariantCaller"))
			{
				try
				{
					DefaultVariantCaller =
						(SupportedVariantCallers)Enum.Parse(typeof(SupportedVariantCallers),
						ParseSettingValue("DefaultVariantCaller"), true);
				}
				catch (Exception)
				{
					// not a valid setting, ignore
				}
			}

			DisablePlotGeneration = GetBooleanSetting("DisablePlotGeneration", DisablePlotGeneration);
			EnableBamStatsByProcess = GetBooleanSetting("EnableBamStatsByProcess", EnableBamStatsByProcess);

			FilterNonPFReads = GetBooleanSetting("FilterNonPFReads", FilterNonPFReads);

			if (_dictionaryConfig.ContainsKey("GATKDownsampleDepth"))
			{
				GATKDownsampleDepth = int.Parse(ParseSettingValue("GATKDownsampleDepth"));
			}

			if (_dictionaryConfig.ContainsKey("GenomePath"))
			{
				GenomePath = ParseSettingValue("GenomePath");
			}

			if (_dictionaryConfig.ContainsKey("IndelRepeatFilterCutoff"))
			{
				IndelRepeatFilterCutoff = int.Parse(ParseSettingValue("IndelRepeatFilterCutoff"));
			}

			if (_dictionaryConfig.ContainsKey("MaximumMemoryGB"))
			{
				MaximumMemoryGB = float.Parse(ParseSettingValue("MaximumMemoryGB"));
			}
			if (_dictionaryConfig.ContainsKey("MaximumHoursPerProcess"))
			{
				MaximumHoursPerProcess = float.Parse(ParseSettingValue("MaximumHoursPerProcess"));
			}
			if (_dictionaryConfig.ContainsKey("MaximumMegabasesIsaacIndex"))
			{
				MaximumMegabasesIsaacIndex = int.Parse(ParseSettingValue("MaximumMegabasesIsaacIndex"));
			}
			if (_dictionaryConfig.ContainsKey("NMaskShortAdapterReads"))
			{
				NMaskShortAdapterReads = int.Parse(ParseSettingValue("NMaskShortAdapterReads"));
			}

			if (_dictionaryConfig.ContainsKey("MinimumAlignReadLength"))
			{
				MinimumAlignReadLength = int.Parse(ParseSettingValue("MinimumAlignReadLength"));
				if (MinimumAlignReadLength < 9 || MinimumAlignReadLength > 50)
				{
					throw new Exception(string.Format("Invalid setting {0} for MinimumAlignReadLength", MinimumAlignReadLength));
				}
			}

			ReportDemuxToInterops = GetBooleanSetting("ReportDemuxToInterops", ReportDemuxToInterops);

			if (_dictionaryConfig.ContainsKey("PostRunCommand"))
			{
				PostRunCommand = ParseSettingValue("PostRunCommand");
			}
			RequirePerfectMiRNAHit = GetBooleanSetting("RequirePerfectMiRNAHit", RequirePerfectMiRNAHit);
			RetainTempFiles = GetBooleanSetting("RetainTempFiles", RetainTempFiles);
			ReportIndelRepeatCount = GetBooleanSetting("ReportIndelRepeatCount", ReportIndelRepeatCount);
			RunGATKForSingleSamples = GetBooleanSetting("RunGATKForSingleSamples", RunGATKForSingleSamples);
            Utilities.RunProcessesWithLowPriority = GetBooleanSetting("RunProcessesWithLowPriority", Utilities.RunProcessesWithLowPriority);

			if (_dictionaryConfig.ContainsKey("SVAnnotation"))
			{
				string tempString = ParseSettingValue("SVAnnotation").ToLowerInvariant();
				switch (tempString)
				{
					case "nirvana":
					case "iae":
					case "illumina annotation engine":
					case "illuminaannotationengine":
						SVAnnotation = SupportedVariantAnnotators.Nirvana;
						break;
					case "svannotate":
						SVAnnotation = SupportedVariantAnnotators.SVAnnotate;
						break;
					case "none":
						SVAnnotation = SupportedVariantAnnotators.None;
						break;
					default:
						throw new Exception("Error: Invalid setting for SVAnnotation. Choose SVAnnotate or Nirvana or None");
				}
			}

			if (_dictionaryConfig.ContainsKey("TempFolder"))
			{
				TempFolder = ParseSettingValue("TempFolder");
			}

			UsePipedBwa = GetBooleanSetting("UsePipedBwa", UsePipedBwa);

			if (_dictionaryConfig.ContainsKey("VariantAnnotation"))
			{
				string tempString = ParseSettingValue("VariantAnnotation").ToLowerInvariant();
				switch (tempString)
				{
					case "nirvana":
					case "iae":
					case "illumina annotation engine":
					case "illuminaannotationengine":
						SVAnnotation = SupportedVariantAnnotators.Nirvana;
						break;
					case "none":
						VariantAnnotation = SupportedVariantAnnotators.None;
						break;
					default:
						throw new Exception("Error: Invalid setting for VariantAnnotation");
				}
			}
			if (_dictionaryConfig.ContainsKey("VariantFilterQualityCutoff"))
			{
				VariantFilterQualityCutoff = int.Parse(ParseSettingValue("VariantFilterQualityCutoff"));
			}
			//this.CheckForUnknownSettings();
		}

		/// <summary>
		///     This constructor is private, to enforce a singleton IsasConfiguration object.
		///     Access it by calling GetConfiguration.
		/// </summary>
		private IsasConfiguration()
			: this(GetAppSettingDictionary())
		{
		}

		/// <summary>
		///     Get the configuration by calling GetConfiguration
		/// </summary>
		public static IsasConfiguration GetConfiguration()
		{
			lock (ConfigLock)
			{
				if (Config == null)
				{
					Config = new IsasConfiguration();
				}
			}
			return Config;
		}

		/// <summary>
		/// Gets the configuration object by providing a Dictionary of the config settings 
		/// </summary>
		/// <param name="settings"></param>
		/// <returns></returns>
		public static IsasConfiguration GetConfiguration(Dictionary<string, string> settings)
		{
			lock (ConfigLock)
			{
				if (Config == null)
				{
					Config = new IsasConfiguration(settings);
				}
			}
			return Config;
		}

		/// <summary>
		/// Returns the app settings as a Dictionary
		/// </summary>
		/// <returns>The app setting laid into a Dictionary</returns>
		public static Dictionary<string, string> GetAppSettingDictionary()
		{
			Dictionary<string, string> result = new Dictionary<string, string>();

			if (!ConfigurationManager.AppSettings.HasKeys()) return result;
			string[] keys = ConfigurationManager.AppSettings.AllKeys;
			foreach (string key in keys)
			{
				string[] values = ConfigurationManager.AppSettings.GetValues(key);
				if (values != null)
				{
					foreach (string value in values)
					{
						result.Add(key, value);
					}
				}
			}
			return result;
		}

		/// <summary>
		///     Get a true/false value from the app.config.  Handle the case where the value isn't set, where it's
		///     set to "1" or "0", or wheere it's set to "true" or "false".  Integer values other than 1 are considered
		///     to be false.
		/// </summary>
		private bool GetBooleanSetting(string name, bool defaultSetting)
		{
			if (!_dictionaryConfig.ContainsKey(name)) return defaultSetting;
			string value = ParseSettingValue(name);
			int intValue;
			bool result = int.TryParse(value, out intValue);
			if (result)
			{
				return intValue == 1;
			}
			bool boolSettingValue;
			result = bool.TryParse(value, out boolSettingValue);
			if (result) return boolSettingValue;

			// The setting is not a valid int or a valid bool.  What is it?  We don't know!  Report an exception!
			throw new Exception(string.Format("Error parsing setting '{0}' = '{1}' from configuration file", name, value));
		}

		/// <summary>
		/// Returns the value associated to a given key in the _dictionaryConfig object
		/// </summary>
		/// <param name="key">Key used to return value</param>
		/// <returns>The value associated with the key</returns>
		private string ParseSettingValue(string key)
		{
			string value;
			_dictionaryConfig.TryGetValue(key, out value);
			if (value != null)
				value = value.Trim();
			if (string.IsNullOrEmpty(value))
			{
				throw new Exception(string.Format("Error: {0} is not properly defined in config settings", key));
			}
			return value;
		}

        public static string GetVersionString()
        {
            return FileVersionInfo.GetVersionInfo(Assembly.GetExecutingAssembly().Location).ProductVersion;
        }
    }
}
