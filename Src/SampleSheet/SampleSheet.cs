// **** IMPORTANT NOTE:
// This source file is used by several projects: Isis, MSR, RTA, HCS, MCS, and BaseSpace
// When making changes, to limit the need to merge, please make fixes/edits first in the Illumina.Isis trunk,
// then propagate the changes to other branches.  
// ****

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Xml.Serialization;
using ProtoBuf;
using Illumina.Common;

namespace SequencingFiles
{
	[ProtoContract]
	public class SampleSheet
	{
		#region serializable members
		// ReSharper disable InconsistentNaming
		[ProtoMember(1)]
		public SampleSheetType Type;
		[ProtoMember(2)]
		public string WorkflowType;
		[ProtoMember(3)]
		public SampleSheetHeader Header { get; set; }
		// ReSharper restore InconsistentNaming
		private ErrorHandler _error;
		private ErrorHandler _log;
		#endregion

		#region members
		// contains manifests if applicable
		[XmlIgnore]
		[ProtoMember(4)]
		public Dictionary<string, string> ManifestLookup = new Dictionary<string, string>();

		[XmlIgnore]
		[ProtoMember(5)]
		public List<string> Manifests
		{
			get { return ManifestLookup.Values.ToList(); }
		}

		[XmlIgnore]
		[ProtoMember(6)]
		public List<Sample> Samples
		{
			get { return _samples; }
		}

		// contains settings strings which depend on workflow.  Keys are lowercase.
		[XmlIgnore]
		[ProtoMember(8)]
		public Dictionary<string, string> Settings { get; set; }

		// List of samples per lane and per index 
		[XmlIgnore]
		[ProtoMember(9)]
		public List<int> ReadCycles { get; set; }

		// Expose the sample sheet sections
		[XmlIgnore]
		[ProtoMember(10)]
		public Dictionary<string, List<string>> Sections { get { return _sections; } }

		// Expose the sample columns
		[XmlIgnore]
		[ProtoMember(11)]
		public Dictionary<string, int> SampleColumns { get { return _sampleColumns; } }

		// raw lines by section in file
		protected readonly Dictionary<string, int> _sampleColumns = new Dictionary<string, int>(StringComparer.InvariantCultureIgnoreCase);
		private readonly List<Sample> _samples = new List<Sample>();
		private bool _requireSample = true;
		private readonly Dictionary<string, List<string>> _sections = new Dictionary<string, List<string>>();

		public static object SampleSheetLock = new object();
		private static readonly char[] HeaderSplitChars = new[] { '[', ']', ',' };
		public static readonly string NoGenomeFolder = "No Path";

		private readonly char[] _slashChars = new[] { '/', '\\' };
		#endregion

		#region header keys

		// these should be in lowercase
		public const string InvestigatorHeader = "Investigator Name";
		public const string ExpNameHeader = "Experiment Name";
		public const string DateHeader = "Date";
		public const string ChemistryHeader = "Chemistry";
		public const string WorkflowHeader = "Workflow";
		public const string ApplicationHeader = "Application";
		public const string AssayHeader = "Assay";

		#endregion

		// constructor
		public SampleSheet()
		{
			Header = new SampleSheetHeader();
			Settings = new Dictionary<string, string>(StringComparer.InvariantCultureIgnoreCase);
			ReadCycles = new List<int>();
		}

		public static SampleSheet DeepClone(SampleSheet sheet)
		{
			return DeepClone(sheet, typeof(SampleSheet.Sample));
		}

		public static SampleSheet DeepClone(SampleSheet sheet, Type sampleType)
		{
			// build protobuf model to take into account potentially extended Sample
			ProtoBuf.Meta.RuntimeTypeModel model = ProtoBuf.Meta.TypeModel.Create();

			var type = model.Add(typeof(SampleSheet), true);
			Type sampleSheetType = sheet.GetType();
			if (sampleSheetType != typeof(SampleSheet))
			{
				// sampleSheetType extended SampleSheet
				model.Add(sampleSheetType, true);
				type.AddSubType(49, sampleSheetType);
			}

			type = model.Add(typeof(SampleSheet.Sample), true);
			if (sampleType != typeof(SampleSheet.Sample))
			{
				// sampleType extended SampleSheet.Sample
				model.Add(sampleType, true);
				type.AddSubType(50, sampleType);
			}

			// serialize and deserialize
			using (MemoryStream stream = new MemoryStream())
			{
				model.Serialize(stream, sheet);
				stream.Position = 0;
				SampleSheet newsheet = null;
				newsheet = (SampleSheet)model.Deserialize(stream, newsheet, typeof(SampleSheet));
				return newsheet;
			}
		}

		public delegate void ErrorHandler(string message);

		/// <summary>
		///     writes to the standard error log
		/// </summary>
		protected void OnError(string message)
		{
			if (_error != null) _error(message);
			else Console.WriteLine(message);
		}

		/// <summary>
		///     writes to the standard output log
		/// </summary>
		protected void OnLog(string message)
		{
			if (_log != null) _log(message);
			else Console.WriteLine(message);
		}

		#region accessory methods for settings

		/// <summary>
		/// returns true if the setting is present
		/// </summary>
		public bool GetBooleanSetting(string settingName, out bool value)
		{
			value = false;
			if (!Settings.ContainsKey(settingName)) return false;

			// Stathis bug fix: Trim flanking space!
			string settingValue = Settings[settingName].ToLowerInvariant().Trim();

			bool isFalse = TrueFalse.IsFalse(settingValue);
			bool isTrue = TrueFalse.IsTrue(settingValue);

			// If it's not set to anything we recognize, report an error!
			if (!isFalse && !isTrue)
			{
				throw new Exception(string.Format("Errror: Sample sheet setting '{0} = '{1}' not understood; use true or false (or 1 or 0) for this setting", settingName, settingValue));
			}

			if (isTrue) value = true;

			// Return true - setting is present
			return true;
		}

		/// <summary>
		///     Returns true if the setting is present.
		/// </summary>
		public bool GetFloatSetting(string settingName, out float value)
		{
			value = 0;
			if (!Settings.ContainsKey(settingName)) return false;

			value = float.Parse(Settings[settingName]);
			return true; // Return true - setting is present
		}

		/// <summary>
		///     Returns true if the setting is present.
		/// </summary>
		public bool GetIntegerSetting(string settingName, out int value)
		{
			value = 0;
			if (!Settings.ContainsKey(settingName)) return false;

			value = int.Parse(Settings[settingName]);
			return true; // Return true - setting is present
		}

		/// <summary>
		/// returns the setting value given a key
		/// </summary>
		public string GetStringSetting(string settingName)
		{
			if (!Settings.ContainsKey(settingName)) return null;
			return Settings[settingName].Trim();
		}

		#endregion

		/// <summary>
		///     Checks the manifest settings in the sample sheet and populates the string builder with
		///     any problems that were discovered during checking.
		/// </summary>
		public void CheckManifestSettings(StringBuilder error, string runFolder)
		{
			// require that the filename without extension is unique even if the files exist in different directories
			// this is so we can write out a unique bed file without name collision
			HashSet<string> fileNames = new HashSet<string>();
			foreach (KeyValuePair<string, string> kvp in ManifestLookup)
			{
				string fileStub = Path.GetFileNameWithoutExtension(kvp.Value);
				if (fileNames.Contains(fileStub))
					error.AppendFormat("Manifest filename must be unique. Found multiple manifest files named {0}.\n", fileStub);
				else
					fileNames.Add(fileStub);
			}

			// make sure we have at least one manifest
			if (ManifestLookup.Count == 0)
			{
				error.AppendFormat("Manifest is required in sample sheet for {0} workflow.\n", WorkflowType);
			}

			// set the manifest filenames for runs with one manifest
			if (ManifestLookup.Count == 1)
			{
				foreach (Sample sample in Samples)
				{
					if (string.IsNullOrEmpty(sample.ManifestFileName)) sample.ManifestFileName = Manifests[0];
				}
			}

			// make sure each sample has a manifest filename
			foreach (Sample sample in Samples)
			{
				if (string.IsNullOrEmpty(sample.ManifestFileName) ||
					(sample.ManifestFileName == Sample.MissingManifestPath))
				{
					error.AppendFormat("Manifest is required for sample: {0}\n", sample.SampleID);
				}
				else if (sample.ManifestFileName.StartsWith(Sample.UnknownManifest))
				{
					error.AppendFormat("Unknown manifest {0} in sample: {1}\n", sample.ManifestFileName.Substring(Sample.UnknownManifest.Length),
						sample.SampleID);
				}
			}
		}

		/// <summary>
		///     Checks the genome paths in the sample sheet 
		///     and throws exception if there are any problems
		/// </summary>
		public void CheckGenomePaths(string genomePath)
		{
			try
			{
				// Confirm that genome paths are correct, handling tricky cases like relative-versus-absolute paths, path separator mixups,
				// or the Chromosomes/WholeGenomeFASTA switcheroo:
				SetGenomes(genomePath);
			}
			catch (Exception ex)
			{
				throw new ApplicationException(string.Format("Error validating Sample Sheet:\n"), ex);
			}
		}

		private string CheckAlternativeGenomePaths(string checkPath, bool AlternativeSlashes)
		{
			string alternativePath = null;

			// If the genome folder is missing, try repairing paths ending in "Chromosomes" or "WholeGenomeFASTA":
			string fileName = Path.GetFileName(checkPath);
			if (fileName == null) return null;

			string folderNameLower = fileName.ToLowerInvariant();
			string directoryName = Path.GetDirectoryName(checkPath);
			if (directoryName == null) return null;

			if (folderNameLower == "chromosomes")
			{
				alternativePath = Path.Combine(directoryName, "WholeGenomeFasta");
				if (!Directory.Exists(alternativePath))
				{
					alternativePath = Path.Combine(directoryName, "WholeGenomeFASTA");
				}
			}
			else if (folderNameLower == "wholegenomefasta")
			{
				alternativePath = Path.Combine(directoryName, "WholeGenomeFasta"); // Try fixing capitalization
				if (!Directory.Exists(alternativePath))
				{
					alternativePath = Path.Combine(directoryName, "WholeGenomeFASTA"); // Try old capitalization
				}
				if (!Directory.Exists(alternativePath))
				{
					alternativePath = Path.Combine(directoryName, "Chromosomes");
				}
			}

			if (Directory.Exists(alternativePath)) return alternativePath;

			if (!AlternativeSlashes)
			{
				// Desperation: Try flip-flopping forward and back-slashes in the path:
				if (checkPath.IndexOf('\\') != -1)
				{
					checkPath = checkPath.Replace('\\', '/');
					if (Directory.Exists(checkPath)) return checkPath;
					alternativePath = CheckAlternativeGenomePaths(checkPath, true);

					if (alternativePath != null) return alternativePath;
				}
				if (checkPath.IndexOf('/') != -1)
				{
					checkPath = checkPath.Replace('/', '\\');
					if (Directory.Exists(checkPath)) return checkPath;
					alternativePath = CheckAlternativeGenomePaths(checkPath, true);
					if (alternativePath != null) return alternativePath;
				}
			}
			return null;
		}

		/// <summary>
		///     Get the array of genome folders (typically just 1) for this sample sheet, given the repository.
		///     (If GenomeRepository is null, then don't complain about directories that don't exist)
		///     This includes some sneaky logic: 
		///     - If we're looking for foo\WholeGenomeFASTA, and it doesn't exist, look for foo\WholeGenomeFasta
		///     - If we're looking for foo\Chromosomes, look for foo\WholeGenomeFasta
		///     - And try flipping back-slashes to forward-slashes (and vice versa)
		///     This sneaky logic is mainly there to smooth how legacy data is handled.  In the past, our references used 
		///     the WholeGenomeFASTA capitalization (and on Linux, this difference matters).  And in the distant past,
		///     our references used the Chromosomes path. 
		/// </summary>
		public void SetGenomes(string genomeRepository)
		{
			List<string> genomes = Sample.GetGenomes(Samples);
			for (int genomeNumber = 1; genomeNumber <= genomes.Count; genomeNumber++)
			{
				foreach (Sample mySample in Samples)
				{
					if (mySample.GenomePath != genomes[genomeNumber - 1]) continue; // Skip over samples for different genomes
					// reset it (relative path will become absolute path).
					// so that developers can call Sheet.Samples[0].GenomePath without worrying about it being relative path
					mySample.GenomePath = GetProperGenomeFolderPath(mySample.GenomePath, genomeRepository);
					mySample.GenomeNumber = genomeNumber;
				}
			}
		}

		/// <summary>
		/// Attach isisConfigGenomeRepository to front of sampleGenomeFolder if the latter is relative path.
		///     This function is split from GetGenomes() so that it could be invoked standalone for one sample.
		/// </summary>
		public string GetProperGenomeFolderPath(string sampleGenomeFolder, string isisConfigGenomeRepository = null)
		{
			string properGenomeFolderPath = null;
			string alternativePath;

			if (sampleGenomeFolder == NoGenomeFolder)
			{
				// "No Path" is a dummy path:
				properGenomeFolderPath = sampleGenomeFolder;
				return properGenomeFolderPath;
			}

			// First, check path relative to repository (if repository is set)
			if (isisConfigGenomeRepository != null)
			{
				string checkPath = Path.Combine(isisConfigGenomeRepository, sampleGenomeFolder);
				if (Directory.Exists(checkPath))
				{
					properGenomeFolderPath = checkPath;
				}
				else
				{
					alternativePath = CheckAlternativeGenomePaths(checkPath, false);
					if (alternativePath != null)
					{
						Console.WriteLine("Switch sample to alt. path {0}", alternativePath);
						properGenomeFolderPath = alternativePath;
					}
				}
			}
			if (properGenomeFolderPath == null)
			{
				// We didn't find this folder below the genome repository.  Assume that it's an absolute path:
				if (Directory.Exists(sampleGenomeFolder))
				{
					properGenomeFolderPath = sampleGenomeFolder;
				}
				else
				{
					// Ok, it doesn't exist as an absolute path either.  Try tweaking the last folder name:
					alternativePath = CheckAlternativeGenomePaths(sampleGenomeFolder, false);
					if (alternativePath != null)
					{
						Console.WriteLine("Switch sample to alt. path {0}", alternativePath);
						properGenomeFolderPath = alternativePath;
					}
				}

			}
			// If we reach this point, we haven't seen a folder at all.  That's not good!
			if (properGenomeFolderPath == null && isisConfigGenomeRepository != null)
			{
				throw new Exception(
					string.Format(
						"Error - GenomeFolder '{0}' was not found, as an absolute path or relative to genome repository '{1}' Please check the sample sheet.",
						sampleGenomeFolder, isisConfigGenomeRepository));
			}
			return properGenomeFolderPath;
		}

		/// <summary>
		/// Parse a sample sheet (either Casava-style or MiSeq-style)
		/// </summary>
		public void ReadSampleSheet(StreamReader reader, string filePath, string runFolder)
		{
			ReadSampleSheet(reader, filePath, runFolder, typeof(Sample));
		}

		public void ReadSampleSheet(StreamReader reader, string filePath, string runFolder, Type sampleType)
		{
			lock (SampleSheetLock)
			{
				bool parsed = false;

				using (reader)
				{
					while (!reader.EndOfStream)
					{
						string line = reader.ReadLine();
						if (string.IsNullOrEmpty(line) || line[0].Equals('#')) continue; // skip blank line or comment

						String[] fields = CSVReader.ParseCommaDelimitedLine(line);

						if (fields[0].ToUpper().Equals("FCID"))
						{
							Type = SampleSheetType.CASAVA;
							ReadCasavaSampleSheet(reader, runFolder);
							parsed = true;
						}
						else
						{
							Type = SampleSheetType.MiSeq;
							ReadMiSeqSampleSheet(reader, line, filePath);
							ParseMiSeqSampleSheet(runFolder, sampleType);
							parsed = true;
						}
					}
				}

				if (!parsed)
				{
					throw new ApplicationException(string.Format("Empty sample sheet at {0}", filePath));
				}

				FinalizeSampleSheet();
			}
		}

		/// <summary>
		///     Subclasses can override this method for additional logic after the entire sample sheet has been read
		/// </summary>
		protected virtual void FinalizeSampleSheet() { }

		/// <summary>
		///     Parse a sample sheet from a file
		/// </summary>
		public void ReadSampleSheet(string filePath, string runFolder)
		{
			ReadSampleSheet(filePath, runFolder, typeof(Sample));
		}
		public void ReadSampleSheet(string filePath, string runFolder, Type sampleType)
		{
			if (!File.Exists(filePath))
			{
				throw new ApplicationException(string.Format("Sample sheet not found at '{0}'", filePath));
			}

			using (FileStream sheetStream = new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.Read))
			using (StreamReader reader = new StreamReader(sheetStream))
			{
				ReadSampleSheet(reader, filePath, runFolder, sampleType);
			}
		}

		/// <summary>
		/// parses a CASAVA sample sheet
		/// </summary>
		private void ReadCasavaSampleSheet(StreamReader reader, string runFolder)
		{
			while (!reader.EndOfStream)
			{
				string line = reader.ReadLine();
				if (string.IsNullOrEmpty(line) || line[0].Equals('#')) continue;

				string[] lineTokens = CSVReader.ParseCommaDelimitedLine(line);

				if (lineTokens != null && lineTokens.Length >= 10)
				{
					int lane = int.Parse(lineTokens[1]);
					string sampleID = lineTokens[2];
					string sampleRef = lineTokens[3];
					string sampleProject = lineTokens[9];
					string index = lineTokens[4];
					string indexSeq1 = "";
					string indexSeq2 = "";

					if (!string.IsNullOrEmpty(index))
					{
						// index 1 and index 2 are separate by '-'
						string[] indexes = index.Split('-');
						indexSeq1 = indexes[0];
						indexSeq2 = "";
						if (indexes.Length > 1) indexSeq2 = indexes[1];
					}

					Sample dna = new Sample
					{
						SampleID = sampleID,
						Name = sampleID,
						Index = _samples.Count + 1,
						SampleProject = sampleProject,
						GenomePath = sampleRef
					};

					dna.RunFolders.Add(runFolder);
					string key = string.Format("{0}|{1}", runFolder, lane);
					dna.RunFoldersAndLanes[key] = string.Format("{0}|{1}", indexSeq1, indexSeq2);
					_samples.Add(dna);
				}
			}

			// that's the only option for legacy CASAVA-format sample sheets
			WorkflowType = "GenerateFASTQ";
		}

		/// <summary>
		/// parses an Isis (non-Casava) sample sheet
		/// </summary>
		private void ReadMiSeqSampleSheet(StreamReader reader, string line, string sampleSheetFilePath)
		{
			Regex SectionHeaderRegex = new Regex(@"^\[(.*)\]");
			while (!reader.EndOfStream)
			{
				bool done = false;
				Match SectionMatch = null;
				if (line != null) SectionMatch = SectionHeaderRegex.Match(line);
				if (SectionMatch != null && SectionMatch.Success)
				{
					string sectionName = SectionMatch.Groups[1].Value;

					if (_sections.ContainsKey(sectionName))
						throw new ApplicationException(string.Format("Sample Sheet: Redundant Section Label in {0}", sampleSheetFilePath));

					_sections.Add(sectionName, new List<string>());

					do
					{
						line = reader.ReadLine();

						// found next header
						SectionMatch = null;
						if (line != null) SectionMatch = SectionHeaderRegex.Match(line);
						if (SectionMatch != null && SectionMatch.Success) done = true;
						else _sections[sectionName].Add(line);
					} while (!done && !reader.EndOfStream);
				}
				else line = reader.ReadLine();
			}
		}

		/// <summary>
		/// Get the full path to the manifest file.  Normally we just look for the
		/// specified file name in the run folder.  Alternatives we support:
		/// - Specified file name lacks .txt extension
		/// - Specified file name lacks .AmpliconManifest extension
		/// - Specified file name is an absolute path (not just a filename)
		/// </summary>
		public static string GetManifestPath(string fileName, string runFolder)
		{
			string finalPath = null;

			List<string> paths = new List<string>();
			paths.Add(Path.Combine(runFolder, Path.GetFileName(fileName)));
			paths.Add(Path.Combine(runFolder, Path.GetFileName(fileName) + ".txt"));
			paths.Add(Path.Combine(runFolder, Path.GetFileName(fileName) + ".ampliconmanifest"));
			paths.Add(Path.Combine(runFolder, Path.GetFileName(fileName) + ".AmpliconManifest"));

			// If it looks like an absolute file path, try treating it like one:
			if (fileName.Contains('/') || fileName.Contains('\\'))
				paths.Add(fileName);

			bool found = false;
			foreach (string path in paths)
			{
				if (File.Exists(path))
				{
					if (!string.IsNullOrEmpty(finalPath) && string.Compare(path, finalPath, true) != 0)
						throw new ArgumentException(string.Format("Ambiguous manifest path defined. Both {0} and {1} exist.\n", finalPath, path));
					finalPath = path;
					found = true;
				}
			}

			if (!found)
			{
				//keep it unchanged, maybe the worker knows the correct directory where the file lives
				return fileName;
			}

			return finalPath;
		}

		/// <summary>
		///     Helper function for ParseMiSeqSampleSheet: Handle [Header] section
		/// </summary>
		private bool ParseSampleSheetHeader()
		{
			if (!_sections.ContainsKey("Header") || _sections["Header"].Count <= 1) return false;
			foreach (string line in _sections["Header"])
			{
				string[] lineTokens = CSVReader.ParseCommaDelimitedLine(line);

				if (lineTokens != null && lineTokens.Length > 1)
				{
					if (lineTokens[0].Equals(InvestigatorHeader, StringComparison.InvariantCultureIgnoreCase))
						Header.InvestigatorName = lineTokens[1];
					else if (lineTokens[0].Equals(ExpNameHeader, StringComparison.InvariantCultureIgnoreCase))
						Header.ExperimentName = lineTokens[1];
					else if (lineTokens[0].Equals(DateHeader, StringComparison.InvariantCultureIgnoreCase))
						Header.Date = lineTokens[1];
					else if (lineTokens[0].Equals(WorkflowHeader, StringComparison.InvariantCultureIgnoreCase))
						WorkflowType = lineTokens[1].Trim(); // trim leading/trailing whitespace from workflow name
					else if (lineTokens[0].Equals(ChemistryHeader, StringComparison.InvariantCultureIgnoreCase))
						Header.Chemistry = lineTokens[1];
					else if (lineTokens[0].Equals(ApplicationHeader, StringComparison.InvariantCultureIgnoreCase))
						Header.Application = lineTokens[1].Trim(); // trim leading/trailing whitespace from application name
					else if (lineTokens[0].Equals(AssayHeader, StringComparison.InvariantCultureIgnoreCase))
						Header.Assay = lineTokens[1].Trim(); // trim leading/trailing whitespace from assay name

				}
			}

			return true;
		}

		/// <summary>
		///     Helper function for ParseMiSeqSampleSheet: Handle [Data] section
		/// </summary>
		private bool ParseDataSection(string runFolder, Type sampleType)
		{
			char[] splitLaneChars = new[] { '+', ';' };

			if (!_sections.ContainsKey("Data") || _sections["Data"].Count <= 1) return false;

			string columns = _sections["Data"][0];
			string[] sampleColumnNames = CSVReader.ParseCommaDelimitedLine(columns);
			int columnIndex = 0;

			foreach (string columnName in sampleColumnNames)
			{
				if (!_sampleColumns.ContainsKey(columnName))
				{
					_sampleColumns.Add(columnName, columnIndex);
				}

				columnIndex++;
			}

			if (!_sampleColumns.ContainsKey("sampleid") && !_sampleColumns.ContainsKey("sample_id"))
			{
				throw new Exception("Sample sheet does not include SampleID column");
			}

			Dictionary<string, Sample> samplesByID = new Dictionary<string, Sample>();

			// Important note: One-based numbering!  Zero is reserved for reads whose index
			// we couldn't resolve.
			int nextSampleIndex = 1;
			const string sampleIdNameCharPattern = @"[^A-Za-z0-9_-]"; // Match anything that's not a letter, number, dash or underscore.  (Also match weird unicode characters that are still considered "letters")
			Regex sampleIdNameRegex = new Regex(sampleIdNameCharPattern);
			Dictionary<string, int> genomeNumbers = new Dictionary<string, int>();

			for (int lineIndex = 1; lineIndex < _sections["Data"].Count; lineIndex++)
			{
				string line = _sections["Data"][lineIndex];
				string[] lineTokens = CSVReader.ParseCommaDelimitedLine(line);

				// Discard empty lines:
				bool blankLine = lineTokens.All(token => token.Trim().Length <= 0);
				if (blankLine) continue;

				// If rightmost columns are missing, add extra lineTokens to the correct length:
				if (lineTokens.Length < _sampleColumns.Count)
				{
					string[] oldTokens = lineTokens;
					lineTokens = new string[_sampleColumns.Count];
					Array.Copy(oldTokens, lineTokens, oldTokens.Length);
				}

				if (lineTokens == null || lineTokens.Length < _sampleColumns.Count) continue;

				// Add a new sample, or look up an existing sample:
				string sampleID = null;
				for (int cIndex = 0; cIndex < _sampleColumns.Count; cIndex++)
				{
					string trimmedString = lineTokens[cIndex] == null ? "" : lineTokens[cIndex].Trim();

					switch (sampleColumnNames[cIndex].ToLowerInvariant())
					{
						case "sample_id":
						case "sampleid":
							if (sampleIdNameRegex.IsMatch(trimmedString))
							{
								throw new ApplicationException(string.Format("Illegal characters in sample ID '{0}' - the ID must consist of letters, numbers, dashes and underscores.", trimmedString));
							}
							sampleID = trimmedString;
							if (string.IsNullOrEmpty(trimmedString))
							{
								throw new ApplicationException("Sample has empty sample ID");
							}
							break;
						case "sample_name":
						case "samplename":
							if (sampleIdNameRegex.IsMatch(trimmedString))
							{
								throw new ApplicationException(string.Format("Illegal characters in sample name '{0}' - the name must consist of letters, numbers, dashes and underscores.", trimmedString));
							}
							break;
					}
				}

				Sample sample;

				if (samplesByID.ContainsKey(sampleID))
				{
					sample = samplesByID[sampleID];
				}
				else
				{
					sample = (Sample)Activator.CreateInstance(sampleType);
					sample.SampleID = sampleID;
					//sample = new Sample { SampleID = sampleID };
					samplesByID[sampleID] = sample;
					sample.Index = nextSampleIndex;
					nextSampleIndex++;
					_samples.Add(sample);
				}

				// By default, assume the DNA's run folder is the same as the run folder where
				// the sample sheet is located:
				string sampleRunFolder = runFolder;
				List<int> tempLanes = new List<int>();

				const string indexCharPattern = "[^ACGTacgtNn]";
				Regex badIndexRegex = new Regex(indexCharPattern);

				string indexSequence1 = "";
				string indexSequence2 = "";

				for (int cIndex = 0; cIndex < _sampleColumns.Count; cIndex++)
				{
					string trimmedString = lineTokens[cIndex] == null ? "" : lineTokens[cIndex].Trim();

					string columnName = sampleColumnNames[cIndex].ToLowerInvariant();
					switch (columnName)
					{
						case "runfolder":
							sampleRunFolder = trimmedString;
							// Remove trailing slashes:
							while (sampleRunFolder.Length > 1 && (sampleRunFolder[sampleRunFolder.Length - 1] == '/' || sampleRunFolder[sampleRunFolder.Length - 1] == '\\'))
							{
								sampleRunFolder = sampleRunFolder.Substring(0, sampleRunFolder.Length - 1);
							}
							break;

						case "sampleid":
						case "sample_id":
							break; // These were already handled

						case "index":
							indexSequence1 = trimmedString;
							if (badIndexRegex.IsMatch(indexSequence1))
								throw new ApplicationException(
									string.Format("Unauthorized characters in index sequence: {0}", trimmedString));
							break;

						case "index2":
							indexSequence2 = trimmedString;
							if (badIndexRegex.IsMatch(indexSequence2))
								throw new ApplicationException(
									string.Format("Unauthorized characters in index sequence: {0}", trimmedString));
							break;

						case "genomefolder":
						case "genome":
							sample.GenomePath = trimmedString.TrimEnd(_slashChars);
							break;

						case "samplename":
						case "sample_name":
							if (sampleIdNameRegex.IsMatch(trimmedString))
								throw new ApplicationException(
									string.Format("Unauthorized characters in Sample Name: {0}", trimmedString));
							sample.Name = trimmedString;
							break;

						case "manifest":
							string manifestFilename;
							if (ManifestLookup.TryGetValue(trimmedString, out manifestFilename))
							{
								sample.ManifestFileName = manifestFilename;
							}
							else
							{
								sample.ManifestFileName = string.IsNullOrEmpty(trimmedString)
									? Sample.MissingManifestPath : string.Format("{0}{1}", Sample.UnknownManifest, trimmedString);
							}
							break;

						case "lane":
						case "lanes":
							string[] lanes = trimmedString.Split(splitLaneChars, StringSplitOptions.RemoveEmptyEntries);

							foreach (string lane in lanes)
							{
								int laneNum;
								try
								{
									laneNum = Convert.ToInt32(lane);
								}
								catch
								{
									throw new ApplicationException(string.Format(
											"Do not understand lane number {0} in sample sheet. Please resolve.", lane));
								}
								if (tempLanes.Contains(laneNum))
								{
									throw new ApplicationException(
										string.Format("Lane number {0} is repeated in sample sheet for the same sample.  Please resolve.", lane));
								}
								tempLanes.Add(laneNum);
							}
							break;

						case "sample_project":
						case "sampleproject":
						case "project":
							sample.SampleProject = trimmedString;
							break;

						case "genesfolder":
							sample.GenesFolder = trimmedString.TrimEnd(_slashChars);
							break;
						case "fastqfolder":
							sample.FastqFolder = trimmedString.TrimEnd(_slashChars);
							break;
						default:
							HandleCustomDataColumn(columnName, trimmedString, sample, lineTokens);
							break;
					}
				}

				// Ensure the sample has a name.  IF no name provided, use the sample ID.:
				if (string.IsNullOrEmpty(sample.Name))
				{
					sample.Name = sample.SampleID;
				}

				// Note the run folder, and lane(s), that this sample is found in:
				if (!sample.RunFolders.Contains(sampleRunFolder))
				{
					sample.RunFolders.Add(sampleRunFolder);
				}

				foreach (int laneNumber in tempLanes)
				{
					string key = string.Format("{0}|{1}", sampleRunFolder, laneNumber);
					if (sample.RunFoldersAndLanes.ContainsKey(key))
					{
						throw new Exception(string.Format("Error: Redundant sample {0} for {1}", sample.SampleID, sampleRunFolder));
					}
					sample.RunFoldersAndLanes[key] = string.Format("{0}|{1}", indexSequence1, indexSequence2);
				}

				if (tempLanes.Count == 0)
				{
					//if there are no lanes specified assume they want them all!
					List<int> lanes = GetAllLanes();
					foreach (int lane in lanes)
					{
						string key = string.Format("{0}|{1}", sampleRunFolder, lane);
						if (sample.RunFoldersAndLanes.ContainsKey(key))
						{
							throw new Exception(string.Format("Error: Redundant sample {0} for {1}", sample.SampleID, sampleRunFolder));
						}
						sample.RunFoldersAndLanes[key] = string.Format("{0}|{1}", indexSequence1, indexSequence2);
					}
				}

				if (string.IsNullOrEmpty(sample.GenomePath)) sample.GenomePath = NoGenomeFolder;
			} // Loop over lines in the [Data] section

			return true;
		}

		/// <summary>
		///    Subclasses can override this if they have custom columns they want to parse
		/// </summary>
		protected virtual void HandleCustomDataColumn(string columnName, string columnValue, Sample sample, string[] lineTokens)
		{
			sample.Data[columnName] = columnValue;
		}

		/// <summary>
		///    Subclasses can override this if they have custom columns they want to write out
		/// </summary>
		public virtual string GetCustomDataColumn(string columnName, Sample sample)
		{
			if (!sample.Data.ContainsKey(columnName)) return null;
			return sample.Data[columnName].Trim();
		}

		private void WriteHeader(StreamWriter writer)
		{
			writer.WriteLine("[Header]");
			if (!string.IsNullOrEmpty(Header.InvestigatorName))
				writer.WriteLine(CSVWriter.GetLine(new string[] { InvestigatorHeader, Header.InvestigatorName }));
			if (!string.IsNullOrEmpty(Header.ExperimentName))
				writer.WriteLine(CSVWriter.GetLine(new string[] { ExpNameHeader, Header.ExperimentName }));
			if (!string.IsNullOrEmpty(Header.Date))
				writer.WriteLine(CSVWriter.GetLine(new string[] { DateHeader, Header.Date }));
			if (!string.IsNullOrEmpty(Header.Chemistry))
				writer.WriteLine(CSVWriter.GetLine(new string[] { ChemistryHeader, Header.Chemistry }));
			if (!string.IsNullOrEmpty(WorkflowType))
				writer.WriteLine(CSVWriter.GetLine(new string[] { WorkflowHeader, WorkflowType }));
			if (!string.IsNullOrEmpty(Header.Application))
				writer.WriteLine(CSVWriter.GetLine(new string[] { ApplicationHeader, Header.Application }));
			if (!string.IsNullOrEmpty(Header.Assay))
				writer.WriteLine(CSVWriter.GetLine(new string[] { ChemistryHeader, Header.Assay }));

		}

		private void WriteReads(StreamWriter writer)
		{
			writer.WriteLine("[Reads]");
			foreach (int cycles in ReadCycles)
			{
				writer.WriteLine(cycles);
			}
		}

		/// <summary>
		///    Helper function to keep the order of the settings from the original sample sheet when we write out the modified one
		/// </summary>
		private int GetSettingKeyIndex(string key)
		{
			if (!_sections.ContainsKey("Settings") || _sections["Settings"] == null) return 0;
			for (int index = 0; index < _sections["Settings"].Count; index++)
			{
				string[] line = CSVReader.ParseCommaDelimitedLine(_sections["Settings"][index]);
				if (line.Length > 0 && line[0].Equals(key, StringComparison.InvariantCultureIgnoreCase))
					return index;
			}
			//if the key is not in the sections it is newly added and should appear at the end
			return int.MaxValue;
		}

		/// <summary>
		///    Helper function to keep the order of the manifests from the original sample sheet when we write out the modified one
		/// </summary>
		private int GetManifestKeyIndex(string key)
		{
			if (!_sections.ContainsKey("Manifests")) return 0;
			for (int index = 0; index < _sections["Manifests"].Count; index++)
			{
				string[] line = CSVReader.ParseCommaDelimitedLine(_sections["Manifests"][index]);
				if (line.Length > 0 && line[0].Equals(key))
					return index;
			}
			//if the key is not in the sections it is newly added and should appear at the end
			return int.MaxValue;
		}

		private void WriteManifests(StreamWriter writer)
		{
			writer.WriteLine("[Manifests]");
			List<string> keyList = new List<string>(ManifestLookup.Keys);
			keyList.Sort((x, y) => GetManifestKeyIndex(x).CompareTo(GetManifestKeyIndex(y)));

			foreach (string key in keyList)
			{
				writer.WriteLine(CSVWriter.GetLine(new string[] { key, ManifestLookup[key] }));
			}
		}

		private void WriteSettings(StreamWriter writer)
		{
			writer.WriteLine("[Settings]");
			List<string> keyList = new List<string>(Settings.Keys);
			keyList.Sort((x, y) => GetSettingKeyIndex(x).CompareTo(GetSettingKeyIndex(y)));

			foreach (string key in keyList)
			{
				writer.WriteLine(CSVWriter.GetLine(new string[] { key, Settings[key] }));
			}
		}

		// if lanes are not specified we assume they want all possible lanes. Currently we assume this means lanes 1-8
		public static List<int> GetAllLanes()
		{
			return Enumerable.Range(1, 8).ToList();
		}

		private void WriteData(StreamWriter writer)
		{
			writer.WriteLine("[Data]");

			List<string> keyList = new List<string>(_sampleColumns.Keys);
			keyList.Sort((x, y) => _sampleColumns[x].CompareTo(_sampleColumns[y]));

			//see if we need Lanes column
			bool requireLanes = false;
			List<int> allLanes = GetAllLanes();
			foreach (Sample s in _samples)
			{
				foreach (string runFolder in s.RunFolders)
				{
					HashSet<int> lanes = new HashSet<int>(s.GetLanesForRunFolder(runFolder));
					if (!lanes.SetEquals(allLanes))
					{
						requireLanes = true;
						break;
					}
				}
				if (requireLanes) break;
			}

			//see if we need Index or Index2 columns            
			bool requireIndex = false;
			bool requireIndex2 = false;
			foreach (Sample s in _samples)
			{
				foreach (KeyValuePair<string, string> kvp in s.RunFoldersAndLanes)
				{
					string[] indexes = kvp.Value.Split('|');
					if (indexes[0] != "") requireIndex = true;
					if (indexes[1] != "") requireIndex2 = true;
					if (requireIndex && requireIndex2) break;
				}
				if (requireIndex && requireIndex2) break;
			}

			//add extra columns if necessary
			if (requireIndex && !keyList.Contains("index", StringComparer.InvariantCultureIgnoreCase))
				keyList.Add("Index");
			if (requireIndex2 && !keyList.Contains("index2", StringComparer.InvariantCultureIgnoreCase))
				keyList.Add("Index2");
			if (requireLanes && !keyList.Contains("lane", StringComparer.InvariantCultureIgnoreCase) && !keyList.Contains("lanes", StringComparer.InvariantCultureIgnoreCase))
				keyList.Add("Lanes");

			//force runfolder column so that there is no ambiguity between the new runfolder where this sample sheet will go and the sample's actual runfolder 
			if (!keyList.Contains("RunFolder", StringComparer.InvariantCultureIgnoreCase))
				keyList.Add("RunFolder");

			//write out the data section column headers
			writer.WriteLine(CSVWriter.GetLine(keyList.ToArray()));

			foreach (Sample s in _samples)
			{
				foreach (string line in GetSampleLines(s, keyList))
					writer.WriteLine(line);
			}
		}

		private IEnumerable<string> GetSampleLines(Sample s, List<string> keyList)
		{
			var lines = new List<string>();

			foreach (string runFolder in s.RunFolders)
			{
				List<string> indexPairs = s.GetIndexPairsForRunFolder(runFolder);
				foreach (string indexPair in indexPairs)
				{
					List<int> lanes = s.GetLanesForRunFolderAndIndex(runFolder, indexPair);
					lanes.Sort();
					string lanesString = "";
					HashSet<int> allLanes = new HashSet<int>(GetAllLanes());
					if (!allLanes.SetEquals(lanes))
						lanesString = string.Join("+", lanes);

					List<string> line = new List<string>();
					foreach (string key in keyList)
					{
						switch (key.ToLowerInvariant())
						{
							case "runfolder":
								line.Add(runFolder);
								break;
							case "sampleid":
							case "sample_id":
								line.Add(s.SampleID);
								break;
							case "index":
								string[] splat = indexPair.Split('|');
								line.Add(splat[0]);
								break;
							case "index2":
								string[] splat2 = indexPair.Split('|');
								line.Add(splat2[1]);
								break;
							case "genomefolder":
							case "genome":
								line.Add(s.GenomePath);
								break;
							case "samplename":
							case "sample_name":
								line.Add(s.Name);
								break;
							case "manifest":
								string manifestKey = ManifestLookup.FirstOrDefault(kvp => kvp.Value == s.ManifestFileName).Key;
								if (manifestKey == null && s.ManifestFileName.StartsWith(Sample.UnknownManifest))
									manifestKey = s.ManifestFileName.Replace(Sample.UnknownManifest, "");
								line.Add(manifestKey ?? string.Empty);
								break;
							case "lane":
							case "lanes":
								line.Add(lanesString);
								break;
							case "sample_project":
							case "sampleproject":
							case "project":
								line.Add(s.SampleProject);
								break;
							case "genesfolder":
								line.Add(s.GenesFolder);
								break;
							case "fastqfolder":
								line.Add(s.FastqFolder);
								break;
							default:
								line.Add(GetCustomDataColumn(key.ToLowerInvariant(), s));
								break;
						}
					}
					lines.Add(CSVWriter.GetLine(line.ToArray()));
				}
			}
			return lines;
		}

		/// <summary>
		///    Write out a sample sheet (current format, not legacy Casava-format) using the parsed data contained in this instance
		/// </summary>
		public void WriteIsisSampleSheet(string path)
		{
			using (StreamWriter writer = new StreamWriter(path))
			{
				//[Header]
				WriteHeader(writer);

				//[Manifests]
				WriteManifests(writer);

				//[Reads]
				WriteReads(writer);

				//[Settings]
				WriteSettings(writer);

				//[Data]
				WriteData(writer);
			}
		}

		/// <summary>
		///     Parse sample sheet data, which has been loaded from a file into _Sections.
		///     Throws an exception if sample sheet is invalid.
		/// </summary>
		private void ParseMiSeqSampleSheet(string runFolder, Type sampleType)
		{
			string[] lineTokens;

			// here is where we put the raw file sections into data structures
			bool haveHeader = ParseSampleSheetHeader();

			if (_sections.ContainsKey("Manifests") && _sections["Manifests"].Count > 0)
			{
				foreach (string line in _sections["Manifests"])
				{
					lineTokens = CSVReader.ParseCommaDelimitedLine(line);
					if (lineTokens != null && lineTokens.Length > 1 && lineTokens[0] != null && lineTokens[0].Length > 0)
					{
						if (ManifestLookup.ContainsKey(lineTokens[0]))
						{
							throw new Exception(string.Format("Error - manifest key {0} was duplicated in sample sheet", lineTokens[0]));
						}
						ManifestLookup.Add(lineTokens[0], SampleSheet.GetManifestPath(lineTokens[1], runFolder));
					}
				}
			}

			if (_sections.ContainsKey("Settings") && _sections["Settings"].Count > 0)
			{
				foreach (string line in _sections["Settings"])
				{
					lineTokens = CSVReader.ParseCommaDelimitedLine(line);
					if (lineTokens != null && lineTokens.Length > 1 && lineTokens[0] != null && lineTokens[0].Length > 0)
					{
						// the Settings dictionary uses case insensitive string comparer for keys so we can keep the original key capitalization for aesthetics
						string Key = lineTokens[0];
						if (Settings.ContainsKey(Key))
						{
							throw new Exception(string.Format("Error - sample sheet setting {0} was duplicated in sample sheet.  Setting names can appear at most once in the [Settings] section.  Some settings allow multiple values; see the documentation for details.", lineTokens[0]));
						}
						Settings.Add(Key, lineTokens[1]);

						// Sanity-checking: We should only have two (non-empty) tokens!)
						for (int tokenIndex = 2; tokenIndex < lineTokens.Length; tokenIndex++)
						{
							if (lineTokens[tokenIndex] != null && lineTokens[tokenIndex].Trim().Length > 0)
							{
								throw new Exception(string.Format("Error: [Settings] line contains more than two columns: '{0}'", line));
							}
						}
					}
				}
			}

			if (_sections.ContainsKey("Reads") && _sections["Reads"].Count > 0)
			{
				for (int lineIndex = 0; lineIndex < _sections["Reads"].Count; lineIndex++)
				{
					string text = _sections["Reads"][lineIndex];
					if (string.IsNullOrEmpty(text) == false)
					{
						int result;
						string readsString = CSVReader.ParseCommaDelimitedLine(_sections["Reads"][lineIndex])[0];

						// Ignore empty lines, and ignore "index" lines for back-compatibility:
						if (string.IsNullOrEmpty(readsString) || readsString == "" ||
							String.Compare(readsString, "index", StringComparison.OrdinalIgnoreCase) == 0)
						{
							continue;
						}

						bool isInt = int.TryParse(readsString, out result);
						if (!isInt)
						{
							throw new Exception(string.Format("Error reading Read Cycles value '{0}' in sample sheet.", readsString));
						}
						ReadCycles.Add(result);
					}
				}
			}

			bool haveData = ParseDataSection(runFolder, sampleType);

			if (!haveHeader) throw new ApplicationException("No header section found in sample sheet");
			if (_requireSample && !haveData) throw new ApplicationException("No data section found in sample sheet");
			if (_requireSample && _samples.Count == 0) throw new ApplicationException("No valid sample records found in sample sheet");

			// Sanity check: Sample IDs must be unique!
			Dictionary<string, bool> sampleIDs = new Dictionary<string, bool>();

			foreach (Sample dna in Samples)
			{
				if (sampleIDs.ContainsKey(dna.SampleID))
				{
					throw new ApplicationException(string.Format("Duplicate sample ID {0} in sample sheet", dna.SampleID));
				}

				sampleIDs[dna.SampleID] = true;
			}
		}

		/// <summary>
		///     Wrapper for LoadSampleSheet - try renaming the file if it's locked.
		/// </summary>
		public static SampleSheet Load(string sampleSheetPath)
		{
			return LoadSampleSheetWithRename(sampleSheetPath);
		}

		public static SampleSheet LoadSampleSheetWithRename(string sampleSheetPath, bool requireSample = true, string runFolder = null)
		{
			return LoadSampleSheetWithRename(sampleSheetPath, typeof(Sample), requireSample, null, null, null, runFolder);
		}

		public static SampleSheet LoadSampleSheetWithRename(string sampleSheetPath, Type sampleType, bool requireSample = true, Type sampleSheetType = null, ErrorHandler error = null, ErrorHandler log = null, string runFolder = null)
		{
			if (string.IsNullOrEmpty(runFolder))
				runFolder = Path.GetDirectoryName(sampleSheetPath);
			SampleSheet sheet;
			try
			{
				sheet = Load(sampleSheetPath, runFolder, sampleType, requireSample, sampleSheetType, error, log);
			}
			catch (IOException ex)
			{
				if (!File.Exists(sampleSheetPath)) throw;

				// It's possible that the sample sheet file is locked by excel.  That has a way of happening...
				// Try making a temporary copy.  Use a unique-ified sample name in case two people call 
				// LoadSampleSheetWithRename at the same time:
				try
				{
					string tempCopyPath = null;
					lock (SampleSheetLock)
					{
						tempCopyPath = Path.GetTempFileName(); // Creates the temp file
						File.Delete(tempCopyPath);
						File.Copy(sampleSheetPath, tempCopyPath);
					}
					sheet = Load(tempCopyPath, runFolder, sampleType, requireSample, sampleSheetType, error, log);
					File.Delete(tempCopyPath);
				}
				catch
				{
					throw ex;
				}
			}
			return sheet;
		}

		/// <summary>
		/// Load a sample sheet from the specified file path.
		/// </summary>
		public static SampleSheet Load(string filePath, string runFolder, bool requireSample = true)
		{
			return Load(filePath, runFolder, typeof(Sample), requireSample);
		}

		public static SampleSheet Load(string filePath, string runFolder, Type sampleType, bool requireSample, Type sampleSheetType = null, ErrorHandler error = null, ErrorHandler log = null)
		{
			if (sampleSheetType == null)
				sampleSheetType = typeof(SampleSheet);
			SampleSheet sheet = (SampleSheet)Activator.CreateInstance(sampleSheetType);
			sheet._error = error;
			sheet._log = log;
			sheet._requireSample = requireSample;
			sheet.ReadSampleSheet(filePath, runFolder, sampleType);
			return sheet;
		}

		/// <summary>
		///     Validate sample sheet before saving
		/// </summary>
		public void Validate(string filePath, String genomeFolder = null)
		{
			StringBuilder error = new StringBuilder();

			// Settings validation has been removed from the sample sheet class.  Which settings are valid and invalid
			// is not information that the humble sample sheet object understands; the collection of settings available
			// and the values they can take on is workflow-specific logic.

			// Data
			ValidateDataSection(error, genomeFolder);

			if (ManifestLookup.Count > 0)
			{
				string runFolder = Path.GetDirectoryName(filePath) ?? string.Empty;

				foreach (KeyValuePair<string, string> pair in ManifestLookup)
				{
					// validate if manifest exists
					if (!File.Exists(Path.Combine(runFolder, pair.Value)) &&
						!File.Exists(Path.Combine(runFolder, pair.Value + ".txt")) &&
						!File.Exists(Path.Combine(runFolder, pair.Value + ".ampliconmanifest")))
					{
						error.AppendFormat("Manifest {0} was not found in {1}.", pair.Value, runFolder);
					}
				}
			}

			// throw exception if there is any error
			if (error.Length > 0) throw new ApplicationException(string.Format("Error validating Sample Sheet: {0}", error));
		}

		/// <summary>
		///     Helper function for Validate(): To validate the [Data] section, with the samples.
		/// </summary>
		private void ValidateDataSection(StringBuilder error, String genomeFolder)
		{
			HashSet<string> sampleIDs = new HashSet<string>();
			HashSet<string> badSamples = new HashSet<string>();

			// List all the RunFoldersAndLanes keys:
			HashSet<string> runFoldersAndLanes = new HashSet<string>();

			foreach (Sample dna in Samples)
			{
				foreach (string key in dna.RunFoldersAndLanes.Keys)
				{
					if (!runFoldersAndLanes.Contains(key))
					{
						runFoldersAndLanes.Add(key);
					}
				}
			}

			foreach (string key in runFoldersAndLanes)
			{
				int indexLength = -1;
				int indexLength2 = -1;
				HashSet<string> indexes = new HashSet<string>();

				foreach (Sample dna in Samples)
				{
					if (!dna.RunFoldersAndLanes.ContainsKey(key)) continue;

					// Check index uniqueness:
					string sampleIndex = dna.RunFoldersAndLanes[key];

					if (indexes.Contains(sampleIndex))
					{
						if (!badSamples.Contains(dna.SampleID))
						{
							error.AppendFormat("Error: Duplicate index {0} for {1}", sampleIndex, key);
							badSamples.Add(dna.SampleID);
						}
					}

					indexes.Add(sampleIndex);
					string[] index = sampleIndex.Split('|');

					// Check index length:
					if (indexLength == -1)
					{
						indexLength = index[0].Length;
						indexLength2 = index[1].Length;
					}
					else
					{
						// Check consistent length:
						if (index[0].Length != indexLength)
						{
							if (!badSamples.Contains(dna.SampleID))
							{
								error.AppendFormat("Inconsistent index lengths for {0}", key);
								badSamples.Add(dna.SampleID);
							}
						}

						if (index[1].Length != indexLength2)
						{
							if (!badSamples.Contains(dna.SampleID))
							{
								error.AppendFormat("Inconsistent index2 lengths for {0}", key);
								badSamples.Add(dna.SampleID);
							}
						}
					}
				}
			}

			for (int sampleIndex = 0; sampleIndex < Samples.Count; sampleIndex++)
			{
				Sample dna = Samples[sampleIndex];

				// check for empty sampleID
				if (string.IsNullOrEmpty(dna.SampleID))
					error.AppendFormat("Sample {0}, sample ID is empty.", (sampleIndex + 1)).AppendLine();

				// check for duplicate sampleID
				if (sampleIDs.Contains(dna.SampleID))
					error.AppendFormat("Duplicate sample ID {0} in sample sheet.", dna.SampleID).AppendLine();
				else
					sampleIDs.Add(dna.SampleID);

				// in case the workflow has no manifest
				if (Manifests.Count > 0 && !string.IsNullOrEmpty(dna.ManifestFileName))
				{
					bool foundValue = ManifestLookup.Any(kvp => kvp.Value == dna.ManifestFileName);

					if (!foundValue)
					{
						error.AppendFormat("Sample {0} manifest file name does not exist in the [Manifests] sections.", dna.SampleID).AppendLine();
					}
				}
			}
		}

		public Sample GetSample(string sampleID)
		{
			foreach (Sample sample in Samples)
			{
				if (sample.SampleID == sampleID)
					return sample;
			}
			return null;
		}

		/// <summary>
		///     Helper for Isaac: Save a legacy-format sample sheet for a run folder.
		/// </summary>
		public static void SaveLegacySampleSheet(List<SampleSheet.Sample> samples, string runFolder, string flowCellID, string referenceFolder,
			string outputPath, bool fastQInput, string wholeGenomeFASTA)
		{
			using (StreamWriter sheetWriter = new StreamWriter(outputPath))
			{
				sheetWriter.WriteLine("FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject");

				foreach (Sample dna in samples)
				{
					if (dna.WholeGenomeFastaPath != wholeGenomeFASTA) continue;
					foreach (string key in dna.RunFoldersAndLanes.Keys)
					{
						string[] bits = key.Split('|');
						if (bits[0] != runFolder) continue;

						string theIndex = string.Empty;

						if (!fastQInput)
						{
							string[] indexes = dna.RunFoldersAndLanes[key].Split('|');
							theIndex = indexes[0];
							if (indexes[1].Length > 0)
							{
								theIndex = string.Format("{0}-{1}", indexes[0], indexes[1]);
							}
						}
						string[] line = new string[] { flowCellID, bits[1], dna.SampleID, referenceFolder, theIndex, dna.Name, "N", "", "", "" };
						sheetWriter.WriteLine(CSVWriter.GetLine(line));
					}
				}
			}
		}

		[ProtoContract]
		public class Sample
		{
			// ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializable members that are sensitive to being changed
			#region members
			[ProtoMember(1)]
			public Dictionary<string, string> Data = new Dictionary<string, string>(StringComparer.InvariantCultureIgnoreCase);
			[ProtoMember(2)]
			public Dictionary<string, string> RunFoldersAndLanes { get; set; } // Keys of the form RunFolder|LaneNumber; values of the form index1|index2
			[ProtoMember(3)]
			public List<string> RunFolders { get; set; }
			[ProtoMember(4)]
			public string GenesFolder { get; set; } // annotation directory in iGenomes for whole genome RNA
			[ProtoMember(5)]
			public int GenomeNumber { get; set; } // 1-based
			[ProtoMember(6)]
			public int Index { get; set; } // internal one-based number
			[ProtoMember(7)]
			public string GenomePath { get; set; }
			[ProtoMember(8)]
			public string ManifestFileName { get; set; } // for targeted workflows
			[ProtoMember(9)]
			public string Name { get; set; }
			[ProtoMember(10)]
			public string SampleID { get; set; }
			[ProtoMember(11)]
			public string SampleProject { get; set; }
			[ProtoMember(12)]
			public string WholeGenomeFastaPath { get; set; }
			[ProtoMember(13)]
			public string FastqFolder { get; set; } // Folder to get fastq files from if fastqInput,useExternal or UseExisting is used.
			[ProtoMember(14)]
			public string ManifestName { get; set; } // for targeted workflows to print pretty name in reports
			// Does this sample use staged fastq files? If yes, all files in the <AnalysisDir>/Fastq/<SampleID> folder will be used, regardless of name
			public bool IsStaged = false;
			#endregion

			#region manifest constants
			public const string MissingManifestPath = "(missing)";
			public const string UnknownManifest = "(unknown)";
			#endregion
			// ReSharper restore InconsistentNaming

			// constructor
			public Sample()
			{
				RunFoldersAndLanes = new Dictionary<string, string>();
				RunFolders = new List<string>();
			}

			//copy ctor - create for use by FastQApp
			public Sample(Sample oldSample)
			{
				Data = new Dictionary<string, string>(oldSample.Data);
				RunFoldersAndLanes = new Dictionary<string, string>(oldSample.RunFoldersAndLanes);
				RunFolders = new List<string>(oldSample.RunFolders);
				GenesFolder = oldSample.GenesFolder;
				GenomeNumber = oldSample.GenomeNumber;
				Index = oldSample.Index;
				GenomePath = oldSample.GenomePath;
				ManifestFileName = oldSample.ManifestFileName;
				Name = oldSample.Name;
				SampleID = oldSample.SampleID;
				SampleProject = oldSample.SampleProject;
				WholeGenomeFastaPath = oldSample.WholeGenomeFastaPath;
			}

			public static List<string> GetSourceRunFolders(List<Sample> samples)
			{
				return samples.SelectMany(s => s.RunFolders).Distinct(StringComparer.InvariantCultureIgnoreCase).ToList();
			}

			private static OrderedDictionary<string, List<Sample>> InternalGetSamplesByGenome(List<Sample> samples)
			{
				var samplesByGenome = new OrderedDictionary<string, List<Sample>>(StringComparer.InvariantCultureIgnoreCase);
				foreach (Sample sample in samples)
				{
					if (!samplesByGenome.ContainsKey(sample.GenomePath))
						samplesByGenome[sample.GenomePath] = new List<Sample>();
					samplesByGenome[sample.GenomePath].Add(sample);
				}
				return samplesByGenome;
			}

			public static Dictionary<string, List<Sample>> GetSamplesByGenome(List<Sample> samples)
			{
				return new Dictionary<string, List<Sample>>(InternalGetSamplesByGenome(samples), StringComparer.InvariantCultureIgnoreCase);
			}

			public static List<string> GetGenomes(List<Sample> samples)
			{
				return InternalGetSamplesByGenome(samples).Keys.ToList();
			}

			public List<int> GetLanesForRunFolder(string runFolder)
			{
				List<int> lanes = new List<int>();
				foreach (string key in RunFoldersAndLanes.Keys)
				{
					string[] splat = key.Split('|');
					if (splat[0] == runFolder)
						lanes.Add(int.Parse(splat[1]));
				}
				return lanes;
			}

			public List<string> GetIndexPairsForRunFolder(string runFolder)
			{
				HashSet<string> indexPairs = new HashSet<string>();
				foreach (KeyValuePair<string, string> kvp in RunFoldersAndLanes)
				{
					string[] keySplat = kvp.Key.Split('|');
					if (keySplat[0] == runFolder)
						indexPairs.Add(kvp.Value);
				}
				return indexPairs.ToList();
			}

			public List<int> GetLanesForRunFolderAndIndex(string runFolder, string indexPair)
			{
				List<int> lanes = new List<int>();
				foreach (KeyValuePair<string, string> kvp in RunFoldersAndLanes)
				{
					string[] keySplat = kvp.Key.Split('|');
					if (keySplat[0] == runFolder && kvp.Value == indexPair)
						lanes.Add(int.Parse(keySplat[1]));
				}
				return lanes;
			}
		}

		[ProtoContract]
		public class SampleSheetHeader
		{
			// ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializable members that are sensitive to being changed
			[ProtoMember(1)]
			public string Chemistry;
			[ProtoMember(2)]
			public string Date;
			[ProtoMember(3)]
			public string ExperimentName;
			[ProtoMember(4)]
			public string InvestigatorName;
			[ProtoMember(5)]
			public string Application;
			[ProtoMember(6)]
			public string Assay;
			// ReSharper restore InconsistentNaming
		}
	}

	[ProtoContract]
	public enum SampleSheetType
	{
		// ReSharper disable InconsistentNaming
		[ProtoEnum]
		CASAVA,
		[ProtoEnum]
		MiSeq
		// ReSharper restore InconsistentNaming
	}

	internal static class TrueFalse
	{
		public static string[] TrueStrings = { "true", "t", "yes", "y", "1" };
		public static string[] FalseStrings = { "false", "f", "no", "n", "0" };

		private static readonly HashSet<string> TrueSet;
		private static readonly HashSet<string> FalseSet;

		static TrueFalse()
		{
			TrueSet = new HashSet<string>(TrueStrings);
			FalseSet = new HashSet<string>(FalseStrings);
		}

		public static bool IsTrue(string s)
		{
			return TrueSet.Contains(s.ToLowerInvariant());
		}

		public static bool IsFalse(string s)
		{
			return FalseSet.Contains(s.ToLowerInvariant());
		}
	}
}