using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Text.RegularExpressions;

namespace SequencingFiles
{
	/// <summary>
	///     Contains the basic read information found in the FASTQ files
	/// </summary>
	public struct BoltRead
	{
		public string Bases;
		public string FlowcellID;
		public ushort FragmentAlignmentQuality;
		public string Index;
		public string InstrumentName;
		public bool IsFiltered;
		public string Lane;
		public string Qualities;
		public int ReadNum;
		public string RunID;
		public string Tile;
		public string UnparsedName;
		public string X;
		public string Y;
		public string Header;
		public string UMI;

        /// <summary>
        ///   Parse this.UnparsedName into the appropriate fields.
        /// </summary>
        /// <returns>true if successful</returns>
        public bool ParseReadName()
        {
            string[] nameBits = UnparsedName.Split(':');
            if (nameBits.Length != 7 && nameBits.Length != 8)
                return false;
            
            InstrumentName = nameBits[0];
            RunID = nameBits[1];
            FlowcellID = nameBits[2];
            Lane = nameBits[3];
            Tile = nameBits[4];
            X = nameBits[5];
            Y = nameBits[6];
            UMI = (nameBits.Length > 7) ? nameBits[7] : null;  // UMI is optional

            return true;
        }

        public static string GetUMIFromUnparsedName(string unparsedName)
        {
            string[] nameBits = unparsedName.Split(':');
            return (nameBits.Length == 8) ? nameBits[7] : null;
        }
	}

    // some constants
    public static class Constants
    {
        public const int FastqOffset = 33;
    }

	public abstract class FileCommon : IDisposable
	{
		#region member variables

		private static readonly Regex FilenameFixPattern = new Regex(@"[-]{2,}", RegexOptions.Compiled);
		protected bool IsDisposed;
		protected bool IsOpen;
		protected string FileName;

		#endregion

		/// <summary>
		///     constructor
		/// </summary>
		protected FileCommon()
		{
			IsOpen = false;
			IsDisposed = false;
		}

		public void Dispose()
		{
			Dispose(true);
			GC.SuppressFinalize(this);
		}

		// destructor
		~FileCommon()
		{
			Dispose(false);
		}

		/// <summary>
		///     Closes the file
		/// </summary>
		public abstract void Close();

		// Implement IDisposable
		protected virtual void Dispose(bool disposing)
		{
			lock (this)
			{
				if (!IsDisposed)
				{
					IsDisposed = true;
					Close();
				}
			}
		}

		// Implement IDisposable

		/// <summary>
		///     Returns a list containing all of the paths in the specified directory
		///     that match the supplied regular expression.
		///     N.B. The regular expression is applied to the filename, not the path
		/// </summary>
		public static List<string> GetMatchingPaths(string directory, Regex regex)
		{
			string[] dirPaths = Directory.GetFiles(directory);
			List<string> matchingPaths = new List<string>();

			foreach (string filePath in dirPaths)
			{
				// ignore filenames that do not match the regex
				string filename = Path.GetFileName(filePath);
				Match match = regex.Match(filename);
				if (!match.Success) continue;

				// add the matching path to the list
				matchingPaths.Add(filePath);
			}

			return matchingPaths;
		}

		/// <summary>
		///     Returns a filename-friendly reference name.  
		/// </summary>
		public static string GetChromosomeFilenameFriendlyName(string name)
		{
			if (name == null) return (null);
			name = name.Replace('|', '!').Replace('.', ','); // aesthetics reasons
			string invalidChars = Regex.Escape(new string(Path.GetInvalidFileNameChars()) + ":");
			string invalidReStr = string.Format(@"[{0}]+", invalidChars);
			return Regex.Replace(name, invalidReStr, "-");
		}

		/// <summary>
		///     Returns the path without the extension
		/// </summary>
		public static string GetPathWithoutExtension(string path)
		{
			return Path.Combine(Path.GetDirectoryName(path), Path.GetFileNameWithoutExtension(path));
		}

		public static string StringArrayToListOfParams(string[] strings)
		{
			StringBuilder sb = new StringBuilder("[");
			foreach (string s in strings)
			{
				sb.Append(s + ",");
			}
			sb.Append("]");
			return sb.ToString();
		}

		public static string[] ListOfParamsToStringArray(string param)
		{
			return param.Split(new[] { ',', '[', ']' }, StringSplitOptions.RemoveEmptyEntries);
		}


		/// <summary>
		///     Opens the file
		/// </summary>
		public abstract void Open(string filename);
	}
}