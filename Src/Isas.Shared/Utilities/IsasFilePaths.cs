using System;
using System.IO;
using System.Text.RegularExpressions;

namespace Isas.Shared
{
	/// <summary>
	///     Various hard-coded file names.
	/// </summary>
	public static class IsasFilePaths
	{
		// file names
		public const string SampleSheetFileName = "SampleSheet.csv";
		public const string UsedSampleSheetFileName = "SampleSheetUsed.csv";
		public const string RunInfoFileName = "RunInfo.xml";
		public const string CheckpointFileName = "Checkpoint.txt";
		public const string CompletedJobFileName = "CompletedJobInfo.xml";
		public const string AnalysisErrorFileName = "AnalysisError.txt";
		public const string AnalysisLogFileName = "AnalysisLog.txt";

		// folder names
	    public const string DataFolderName = "Data";
		public const string IntensitiesFolderName = "Intensities";
		public const string BaseCallsFolderName = "BaseCalls";
	    public const string AnalysisFolderName = "Analysis";
		public const string FastqFolderName = "Fastq";

		/// <summary>
		/// Returns the filename for output files specific to a sample
		/// </summary>
		public static string SampleOutputFileName(string sampleName, int sampleNumber, string ext)
		{
			return string.Format("{0}_S{1}.{2}", sampleName, sampleNumber, ext);
		}

		/// <summary>
		/// Returns the file path for the Demux files
		/// </summary>
		public static string DemuxFilePath(string analysisFolder, int lane, int tile, int runFolderNumber)
		{
			return Path.Combine(GetFastqFolder(analysisFolder),
						string.Format("L00{0}", lane),
						string.Format("s_{0}_{1}_{2}.demux", lane, tile, runFolderNumber));
		}

	    public static string GetLoggingFolder(string analysisFolder)
		{
			return Path.Combine(analysisFolder, "Logging");
		}

	    /// <summary>
		/// Root directory for fastq files (actual files are in subdirectories)
		/// </summary>
		public static string GetFastqFolder(string analysisFolder)
		{
			return Path.Combine(analysisFolder, FastqFolderName);
		}

        public static IDirectoryLocation GetFastqFolder(IDirectoryLocation analysisFolder)
        {
            return analysisFolder.CreateSubdirectory(FastqFolderName);
        }

		/// <summary>
		/// Sample-specific directory containing all fastq files of this sample
		/// </summary>
		public static string GetFastqSampleFolder(string analysisFolder, string sampleID)
		{
			return Path.Combine(GetFastqFolder(analysisFolder), sampleID);
		}

		/// <summary>
		/// utility method for transforming the filename extension (e.g. from .genome.vcf.gz to .vcf)
		/// Using String.Replace method will replace all occurences while this method only replaces the occurence at the end
		/// </summary>
		public static string ReplaceFileNameExtension(this string path, string oldSuffix, string newSuffix)
		{
			if (!path.EndsWith(oldSuffix)) return path;
			return path.Substring(0, path.Length - oldSuffix.Length) + newSuffix;
		}

		/// <summary>
		/// Given a file name (e.g. foo_S1.bam or foo_S1.vcf), derive the sample number (e.g. 1).
		/// Deprecated, since in general we should work forward (given the sample sheet, we *know* the file names
		/// that must be present) and not backward (infer sample numbers from what file names are present)
		/// </summary>
		public static int GetSampleNumberFromPath(string bamPath)
		{
			Match theMatch = Regex.Match(bamPath, @"^(.+)_S(\d+)(.+)$");
			if (theMatch == null || !theMatch.Success) throw new Exception(string.Format("Error: Unable to infer sample number from filename {0}", bamPath));
			return int.Parse(theMatch.Groups[2].Value);
		}

	    public static string VariantConcordance(string analysisFolder)
		{
			return Path.Combine(analysisFolder, "VariantConcordance.xml");
		}

		public static string VcfFileName(string analysisFolder, string sampleName, int sampleNumber)
		{
			return Path.Combine(analysisFolder, SampleOutputFileName(sampleName, sampleNumber, "vcf"));
		}

		public static string CarrierReport(string analysisFolder)
		{
			return Path.Combine(analysisFolder, "MiSeqDxCF139VariantAssay.txt");
		}

		public static string CFDiagnosticReport(string analysisFolder)
		{
			return Path.Combine(analysisFolder, "MiSeqDxCFClinicalSequencingAssay.txt");
		}

		public static string DafFolder(string analysisFolder)
		{
			return Path.Combine(analysisFolder, "DataAccessFiles");
		}

		public static string GetConsensusFileName(string analysisFolder, string manifestName)
		{
			manifestName = Path.GetFileNameWithoutExtension(manifestName);
			return Path.Combine(analysisFolder, "Variants", string.Format("s_G1.{0}.consensus.txt", manifestName));
		}

		/// <summary>
		///     Returns the BAM filename corresponding to the sample number
		/// </summary>
		public static string GetBamPathFromSampleNumber(int sampleNumber, string analysisFolder)
		{
			// initialize
			Regex bamRegex = new Regex(@"^.+_S(\d+).bam$", RegexOptions.Compiled | RegexOptions.IgnoreCase);
			string bamFilePath = string.Empty;

			string[] dirFiles = Directory.GetFiles(analysisFolder);

			foreach (string filePath in dirFiles)
			{
				// ignore non-BAM filenames
				string filename = Path.GetFileName(filePath);
				Match bamMatch = bamRegex.Match(filename);
				if (!bamMatch.Success) continue;

				// check the sample number
				int bamSampleNumber = int.Parse(bamMatch.Groups[1].Value);

				if (bamSampleNumber == sampleNumber)
				{
					bamFilePath = filePath;
					break;
				}
			}

			return bamFilePath;
		}

		public static string TargetedRNASeqFileName(string analysisFolder, int manifestNumber)
		{
			return Path.Combine(analysisFolder, string.Format("TargetedRNASeqGeneExpression_M{0}.tsv", manifestNumber));
		}

		public static string GetClassificationPathFromSampleNumber(int sampleNumber, string analysisFolder)
		{
			// initialize
			Regex fileRegex = new Regex(@"^.+_S(\d+).txt.gz$", RegexOptions.Compiled | RegexOptions.IgnoreCase);
			string filePath = string.Empty;
			string[] dirFiles = Directory.GetFiles(analysisFolder);

			foreach (string path in dirFiles)
			{
				// ignore other filenames
				Match fileMatch = fileRegex.Match(Path.GetFileName(path) ?? "");
				if (!fileMatch.Success) continue;
				// check the sample number
				int fileSampleNumber = int.Parse(fileMatch.Groups[1].Value);
				if (fileSampleNumber == sampleNumber)
				{
					filePath = path;
					break;
				}
			}

			return filePath;
		}
	}
}