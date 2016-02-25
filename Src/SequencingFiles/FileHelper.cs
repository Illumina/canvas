using System;
using System.Collections.Generic;
using System.IO;
using Illumina.Common;
using Illumina.Win32;
using Illumina.Zlib;
using System.Security.Cryptography;

namespace SequencingFiles
{
	public static class FileHelper
	{

		public static byte[] GetMD5Checksum(byte[] data)
		{
			MD5CryptoServiceProvider md5Provider = new MD5CryptoServiceProvider();
			return md5Provider.ComputeHash(data);
		}

		public static byte[] GetMD5Checksum(string filepath)
		{
			MD5CryptoServiceProvider md5Provider = new MD5CryptoServiceProvider();
			using (FileStream stream = new FileStream(filepath, FileMode.Open))
			{
				return md5Provider.ComputeHash(stream);
			}
		}

		public static string GetMD5ChecksumString(string filepath)
		{
			return BitConverter.ToString(GetMD5Checksum(filepath)).Replace("-", "").ToLower();
		}

		public static void SaveBCIFile(string bciFilePath, int[] TileNumbers, int[] ClusterCounts)
		{
			using (FileStream myStream = new FileStream(bciFilePath, FileMode.Create, FileAccess.Write, FileShare.Write))
			using (BinaryWriter writer = new BinaryWriter(myStream))
			{
				for (int tileIndex = 0; tileIndex < TileNumbers.Length; tileIndex++)
				{
					writer.Write(TileNumbers[tileIndex]);
					writer.Write(ClusterCounts[tileIndex]);
				}
			}
		}

		public static void LoadBCIFile(string bciFilePath, out int[] TileNumbers, out int[] ClusterCounts)
		{
			List<int> tileNumberList = new List<int>();
			List<int> clusterCountList = new List<int>();
			using (FileStream myStream = new FileStream(bciFilePath, FileMode.Open, FileAccess.Read, FileShare.Read))
			using (BinaryReader reader = new BinaryReader(myStream))
			{
				while (true)
				{
					if (reader.BaseStream.Position >= reader.BaseStream.Length) break;
					//if (reader.PeekChar() < 0) break; // detect EOF on the Stream using BinaryReader -> THIS MAY FAIL!
					int TileNumber = reader.ReadInt32();
					int ClusterCount = reader.ReadInt32();
					if (TileNumber == 0 && ClusterCount == 0) continue; // Ignore any dummy record with TileNumber = ClusterCount = 0
					if (TileNumber < 0 || ClusterCount < 0) throw new Exception(string.Format("Error: Invalid data in .bci file {0}", bciFilePath)); // Sanity-check
					tileNumberList.Add(TileNumber);
					clusterCountList.Add(ClusterCount);
				}
			}
			TileNumbers = tileNumberList.ToArray();
			ClusterCounts = clusterCountList.ToArray();
		}

		/// <summary>
		///     save calls with file sharing
		/// </summary>
		public static void SaveCalls(byte[] calls, string filePath)
		{
			using (FileStream outputStream = new FileStream(filePath, FileMode.Create, FileAccess.Write, FileShare.Read))
			{
				using (BinaryWriter writer = new BinaryWriter(outputStream))
				{
					writer.Write(calls.Length);
					writer.Write(calls);
				}
			}
		}

		/// <summary>
		///     save calls and scores with file sharing
		/// </summary>
		public static void SaveCallsAndScores(byte[] calls, byte[] scores, string filePath)
		{
			using (FileStream outputStream = new FileStream(filePath, FileMode.Create, FileAccess.Write, FileShare.Read))
			{
				using (BinaryWriter writer = new BinaryWriter(outputStream))
				{
					writer.Write(calls.Length);
					for (int clusterIndex = 0; clusterIndex < calls.Length; clusterIndex++)
					{
						Byte callAndScore;
						if (calls[clusterIndex] != 0)
							callAndScore = (byte)((calls[clusterIndex] & 3) + (scores[clusterIndex] << 2));
						else callAndScore = 0;
						writer.Write(callAndScore);
					}
				}
			}
		}

        public static byte[][] LoadCallsAndScores(string filePath, int clusterCount)
		{
			using (FileStream inputStream = new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.Read))
				return LoadCallsAndScores(inputStream, clusterCount);
		}
		public static byte[][] LoadCallsAndScores(string filePath)
		{
			using (FileStream inputStream = new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.Read))
				return LoadCallsAndScores(inputStream);
		}

		// NB: This function doesn't seem to be called by anything, so doesn't check TransparentCompression
		public static byte[][] LoadCompressedCallsAndScores(string filePath, int clusterCount)
		{
			using (ZInputStream zstr = new ZInputStream(filePath))
				return LoadCallsAndScores(zstr, clusterCount);
		}

		// NB: This function doesn't seem to be called by anything, so doesn't check TransparentCompression
		public static byte[][] LoadCompressedCallsAndScores(string filePath)
		{
			using (ZInputStream zstr = new ZInputStream(filePath))
				return LoadCallsAndScores(zstr);
		}

		public static byte[][] LoadBGZippedCallsAndScores(string filePath, ulong startPosition, int clusterCount, out ulong EndPosition)
		{
			byte[][] callsAndScores = null;
			using (BgzipReader reader = new BgzipReader(filePath))
			{
				if (startPosition > 0)
				{
					reader.Seek(startPosition);
				}
				else
				{
					reader.BlockOffset = 4;
				}
				byte[] buffer = reader.ReadBytes(clusterCount, out EndPosition);

				if (callsAndScores == null)
				{
					callsAndScores = new byte[2][];
					callsAndScores[0] = new byte[clusterCount];
					callsAndScores[1] = new byte[clusterCount];
				}
				for (int clusterIndex = 0; clusterIndex < clusterCount; clusterIndex++)
				{
					if (buffer[clusterIndex] != 0) // can't be uncalled, so quality score must be > 0
						callsAndScores[0][clusterIndex] = (byte)((buffer[clusterIndex] & 3) + 8);
					else
						callsAndScores[0][clusterIndex] = 0; // no-call reserved
					callsAndScores[1][clusterIndex] = (byte)(buffer[clusterIndex] >> 2);
				}
			}
			return callsAndScores;
		}

		/// <summary>
		/// Load a big-bcl file:
		/// </summary>
		public static byte[][] LoadCallsAndScores(Stream inputStream, int clusterCount)
		{
			byte[][] callsAndScores = null;
			byte[] buffer = null;
			int tileStart = 0;
			using (BinaryReader reader = new BinaryReader(inputStream))
			{
				while (tileStart < clusterCount)
				{
					int tileClusterCount = reader.ReadInt32();
					if (tileClusterCount < 0)
					{
						throw new ApplicationException("Bad number of clusters in call file " + tileClusterCount.ToString());
					}

					buffer = reader.ReadBytes(tileClusterCount);
					if (buffer.Length != tileClusterCount)
					{
						throw new ApplicationException(
							string.Format("Error: call file has incorrect length (cluster count {0}, read {1} calls)",
										  tileClusterCount, buffer.Length));
					}

					if (callsAndScores == null)
					{
						callsAndScores = new byte[2][];
						callsAndScores[0] = new byte[clusterCount];
						callsAndScores[1] = new byte[clusterCount];
					}
					for (int clusterIndex = 0; clusterIndex < tileClusterCount; clusterIndex++)
					{
						if (buffer[clusterIndex] != 0) // can't be uncalled, so quality score must be > 0
							callsAndScores[0][tileStart + clusterIndex] = (byte)((buffer[clusterIndex] & 3) + 8);
						else
							callsAndScores[0][tileStart + clusterIndex] = 0; // no-call reserved
						callsAndScores[1][tileStart + clusterIndex] = (byte)(buffer[clusterIndex] >> 2);
					}
					tileStart += tileClusterCount;
				}
			}

			return callsAndScores;
		}


		/// <summary>
		///     Load a .bcl file
		/// </summary>
		public static byte[][] LoadCallsAndScores(Stream inputStream)
		{
			byte[][] callsAndScores = null;
			byte[] buffer = null;
			int clusterCount = 0;

			using (BinaryReader reader = new BinaryReader(inputStream))
			{
				clusterCount = reader.ReadInt32();
				if (clusterCount < 0)
				{
					throw new ApplicationException("Bad number of clusters in call file " + clusterCount.ToString());
				}

				buffer = reader.ReadBytes(clusterCount);
				if (buffer.Length != clusterCount)
				{
					throw new ApplicationException(
						string.Format("Error: call file has incorrect length (cluster count {0}, read {1} calls)",
									  clusterCount, buffer.Length));
				}
			}

			if (buffer != null)
			{
				callsAndScores = new byte[2][];
				callsAndScores[0] = new byte[clusterCount];
				callsAndScores[1] = new byte[clusterCount];

				for (int clusterIndex = 0; clusterIndex < clusterCount; clusterIndex++)
				{
					if (buffer[clusterIndex] != 0) // can't be uncalled, so quality score must be > 0
						callsAndScores[0][clusterIndex] = (byte)((buffer[clusterIndex] & 3) + 8);
					else
						callsAndScores[0][clusterIndex] = 0; // no-call reserved
					callsAndScores[1][clusterIndex] = (byte)(buffer[clusterIndex] >> 2);
				}
			}
			return callsAndScores;
		}

		public static void SaveFloatAsText(float[] numbers, string filePath)
		{
			FileInfo fi = new FileInfo(filePath);
			EnsureDirectoryExists(fi.DirectoryName);

			using (StreamWriter writer = new StreamWriter(filePath))
			{
				foreach (float value in numbers)
					writer.WriteLine(value);
			}
		}

		public static float[] LoadTextAsFloat(string filePath)
		{
			List<float> fvals = new List<float>();

			if (File.Exists(filePath))
			{
				string[] fileLines = null;
				using (StreamReader reader = new StreamReader(filePath))
				{
					fileLines = reader.ReadToEnd().Split('\n');
				}
				foreach (string fileLine in fileLines)
				{
					if (fileLine.Length == 0) continue;
					fvals.Add(float.Parse(fileLine));
				}
			}
			return fvals.ToArray();
		}

        //public static void SaveTransformation(Transformation affineTransformation, string filePath)
        //{
        //    if (affineTransformation == null) return;
        //    using (Stream outputStream = new FileStream(filePath, FileMode.Create))
        //    {
        //        using (BinaryWriter writer = new BinaryWriter(outputStream))
        //        {
        //            for (int row = 0; row < 3; row++)
        //            {
        //                for (int column = 0; column < 3; column++)
        //                {
        //                    writer.Write(affineTransformation.Elements[row, column]);
        //                }
        //            }
        //        }
        //    }
        //}

        //public static Transformation LoadTransformation(string filePath)
        //{
        //    if (filePath == null || !File.Exists(filePath))
        //    {
        //        return null;
        //    }
        //    Transformation affineTransformation = new Transformation();
        //    using (Stream inputStream = new FileStream(filePath, FileMode.Open, FileAccess.Read))
        //    {
        //        using (BinaryReader reader = new BinaryReader(inputStream))
        //        {
        //            for (int row = 0; row < 3; row++)
        //            {
        //                for (int column = 0; column < 3; column++)
        //                {
        //                    affineTransformation.Elements[row, column] = reader.ReadDouble();
        //                }
        //            }
        //        }
        //    }
        //    return affineTransformation;
        //}

		public static void SaveAlignmentFile(string filename, short[] alignments)
		{
			int numClusters = alignments.Length;
			using (FileStream fs = new FileStream(filename, FileMode.Create))
			using (BinaryWriter bw = new BinaryWriter(fs))
			{
				bw.Write(numClusters);
				FastWrite(fs, alignments);
			}
		}

		public static short[] LoadAlignmentFile(string filename)
		{
			short[] outputdata = null;
			if (!File.Exists(filename)) return null;
			using (FileStream stream = new FileStream(filename, FileMode.Open))
			using (BinaryReader reader = new BinaryReader(stream))
			{
				int clusterCount = reader.ReadInt32();
				outputdata = FastReadShort(stream, clusterCount);
			}
			return outputdata;
		}

		public static void SaveReadErrorsFile(string filename, short[] readErrors)
		{
			SaveAlignmentFile(filename, readErrors); // same format
		}

		public static short[] LoadReadErrorsFile(string filename)
		{
			return LoadAlignmentFile(filename); // same format
		}

		public static void SaveControlFile(string filename, ushort[] controls)
		{
			int clusterCount = controls.Length;
			const int id1 = 0; // the new controls file have 4 zero bytes, followed by a 4 byte id
			const int id2 = 2;
			using (FileStream stream = new FileStream(filename, FileMode.Create))
			using (BinaryWriter writer = new BinaryWriter(stream))
			{
				writer.Write(id1);
				writer.Write(id2);
				writer.Write(clusterCount);
				FastWrite(stream, controls);
			}
		}

		public static ushort[] LoadControlFile(string filePath)
		{
			ushort[] outputData = null;
			if (!File.Exists(filePath))
				throw new ApplicationException("Control file not found: " + filePath);
			try
			{
				using (FileStream stream = new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.Read))
				using (BinaryReader reader = new BinaryReader(stream))
				{
					int int1 = reader.ReadInt32();
					int int2 = reader.ReadInt32();

					if (int1 == 0 && int2 != 2)
					{
						throw new ApplicationException("Unrecognized filter file id tag.");
					}
					if (int1 != 0)
					{
						// reset reader to beginning, and then parse old format filter files
						reader.BaseStream.Seek(0, SeekOrigin.Begin);
					}
					int clusterCount = reader.ReadInt32();
					outputData = FastReadUshort(stream, clusterCount);
					if (outputData.Length != clusterCount)
						throw new ApplicationException("Incorrect number of clusters in control file" + filePath);
				}
			}
			catch (System.IO.EndOfStreamException)
			{
				throw new ApplicationException("Truncated control file " + filePath);
			}
			return outputData;
		}

		/// <summary>
		///     This will look at a file name and make sure that the directory structure
		///     exist so that a copy can happen. It is expected that the file does not exist
		///     yet
		/// </summary>
		public static bool EnsureFilesDirectoryExists(string filename)
		{
			FileInfo fi = new FileInfo(filename);
			return EnsureDirectoryExists(fi.Directory);
		}

		/// <summary>
		///     This will look at a file name and make sure that the directory structure
		///     exist so that a copy can happen. It is expected that the file does not exist
		///     yet
		/// </summary>
		public static bool EnsureDirectoryExists(string pathName)
		{
			return EnsureDirectoryExists(new DirectoryInfo(pathName));
		}

		/// <summary>
		///     This will look at a directory name and make sure that the directory (and its parents)
		///     exist so that a copy can happen.
		/// </summary>
		public static bool EnsureDirectoryExists(DirectoryInfo di)
		{
			try
			{
				if (di.Exists)
					return true;
				di.Create();
				di = new DirectoryInfo(di.FullName);
				if (di.Exists)
					return true;
			}
			catch (Exception exp)
			{
				Console.WriteLine("Unable to create directory " + exp);
			}
			try
			{
				if (di.FullName.IndexOf(Path.DirectorySeparatorChar) < 0)
					return true; // we are at the root
				EnsureDirectoryExists(di.Parent);
				if (di.Exists)
					return true;
			}
			catch (Exception exp)
			{
				Console.WriteLine("Unable to create parent directory " + exp);
			}
			return false; // unable to create the directory
		}

		private static string GetIntensityFileName(int cycleIndex, int laneIndex, int tileIndex, int channelIndex)
		{
			return string.Format("I_{0}_{1}_{2}_{3}.txt", cycleIndex, laneIndex, tileIndex, channelIndex);
		}

		private static string GetFocusFileName(int cycleIndex, int laneIndex, int tileIndex, int channelIndex)
		{
			return string.Format("F_{0}_{1}_{2}_{3}.txt", cycleIndex, laneIndex, tileIndex, channelIndex);
		}

		private static string GetImageStatsFileName(int cycleIndex, int laneIndex, int tileIndex, int channelIndex)
		{
			return string.Format("C_{0}_{1}_{2}_{3}.txt", cycleIndex, laneIndex, tileIndex, channelIndex);
		}

		public static string[] GetImageStatsFiles(string runFolder)
		{
			string imageStatsDirectory = GetImageStatsDirectory(runFolder);

			if (Directory.Exists(imageStatsDirectory))
				return Directory.GetFiles(imageStatsDirectory, "C_*_*_*_*.txt");

			return new string[0];
		}

		private static int[] ParseImageStatsFileName(string filePath)
		{
			string fileName = (new FileInfo(filePath)).Name;

			int[] indices = new int[4];

			string[] tokens = fileName.Split(new[] { '_', '.' });
			if (tokens.Length != 6) return null; // wrong format

			for (int tokenIndex = 1; tokenIndex < 5; tokenIndex++)
			{
				int index;

				if (!Int32.TryParse(tokens[tokenIndex], out index))
					return null; // bad format

				indices[tokenIndex - 1] = index;
			}

			return indices;
		}

		/// <summary>
		///     Load the focus value from the file saved by SCS
		/// </summary>
		public static float LoadFocusValue(int cycle, int lane, int tile, int channelIndex, string runFolder)
		{
			float focusValue = float.NaN;
			try
			{
				string directory = Path.Combine(Path.Combine(runFolder, "Processed"), "Focus");
				if (!Directory.Exists(directory))
					return focusValue;
				string fileName = GetFocusFileName(cycle, lane, tile, channelIndex);
				string fullFileName = Path.Combine(directory, fileName);
				if (!File.Exists(fullFileName))
					return focusValue;
				using (TextReader reader = new StreamReader(fullFileName))
				{
					focusValue = Convert.ToSingle(reader.ReadLine());
				}
				File.Delete(fullFileName); // we are done with it
			}
			catch (Exception exp)
			{
				Console.WriteLine(exp.ToString());
			}
			return focusValue;
		}

		public static string GetImageStatsDirectory(string runFolder)
		{
			return Path.Combine(Path.Combine(runFolder, "Processed"), "Contrast");
		}

		public static int[] LoadImageStats(string filepath, out int[] tileIndices)
		{
			int[] values = new int[2];
			tileIndices = ParseImageStatsFileName(filepath);

			try
			{
				if (!File.Exists(filepath))
					return null;

				using (TextReader reader = new StreamReader(filepath))
				{
					values[0] = Convert.ToInt32(reader.ReadLine());
					values[1] = Convert.ToInt32(reader.ReadLine());
				}
				File.Delete(filepath); // we are done with it
			}
			catch (Exception exp)
			{
				Console.WriteLine(exp.ToString());
			}
			return values;
		}

		public static void SaveIntensityStats(int cycleIndex, int laneIndex, int tileIndex, int channelIndex,
											  string runFolder, IntensityStats stats)
		{
			string directory = Path.Combine(Path.Combine(runFolder, "Processed"), "Intensity");
			EnsureDirectoryExists(new DirectoryInfo(directory));
			string fileName = GetIntensityFileName(cycleIndex, laneIndex, tileIndex, channelIndex);
			string filePath = Path.Combine(directory, fileName);
			stats.Save(filePath);
		}

		public static IntensityStats LoadIntensityStats(int cycle, int lane, int tile, int channelIndex,
														string runFolder)
		{
			IntensityStats stats = null;
			try
			{
				string directory = Path.Combine(Path.Combine(runFolder, "Processed"), "Intensity");
				if (!Directory.Exists(directory))
					return null;
				string fileName = GetIntensityFileName(cycle, lane, tile, channelIndex);
				string fullFileName = Path.Combine(directory, fileName);
				if (!File.Exists(fullFileName))
					return null;
				stats = IntensityStats.Load(fullFileName);
				File.Delete(fullFileName); // we are done with it
			}
			catch (Exception exp)
			{
				Console.WriteLine(exp.ToString());
			}
			return stats;
		}

		public static void SaveFocusValue(int cycleIndex, int laneIndex, int tileIndex, int channelIndex,
										  string runFolder, float fwhm)
		{
			string focusFilePath = Path.Combine(Path.Combine(runFolder, "Processed"), "Focus");
			EnsureDirectoryExists(new DirectoryInfo(focusFilePath));
			string focusFileName = GetFocusFileName(cycleIndex, laneIndex, tileIndex, channelIndex);
			focusFilePath = Path.Combine(focusFilePath, focusFileName);
			using (TextWriter focusOutput = new StreamWriter(focusFilePath))
			{
				focusOutput.WriteLine(fwhm);
			}
		}

		public static void SaveImageStats(int cycleIndex, int laneIndex, int tileIndex, int channelIndex,
										  string runFolder, int min, int max)
		{
			string statsFilePath = Path.Combine(Path.Combine(runFolder, "Processed"), "Contrast");
			EnsureDirectoryExists(new DirectoryInfo(statsFilePath));
			string statsFileName = GetImageStatsFileName(cycleIndex, laneIndex, tileIndex, channelIndex);
			statsFilePath = Path.Combine(statsFilePath, statsFileName);
			using (TextWriter output = new StreamWriter(statsFilePath))
			{
				output.WriteLine(min);
				output.WriteLine(max);
			}
		}

		/// <summary>
		///     Save passed-filters flags and control flags to a .filter file
		/// </summary>
		public static void SaveFilterFile(string filePath, bool[] passedFilters)
		{
			int clusterCount = passedFilters.Length;
			const int id1 = 0; // the new filter files have 4 zero bytes, followed by a 4 byte id
			const int id2 = 3;
			using (FileStream outputStream = new FileStream(filePath, FileMode.Create))
			using (BinaryWriter writer = new BinaryWriter(outputStream))
			{
				// Filter files have 4 0 bytes, followd by an id int
				writer.Write(id1);
				writer.Write(id2);
				writer.Write(clusterCount);
				for (int clusterIndex = 0; clusterIndex < clusterCount; clusterIndex++)
				{
					byte value = (byte)(passedFilters[clusterIndex] ? 1 : 0);
					writer.Write(value);
				}
			}
		}

		/// <summary>
		///     Load passed-filters flags and control flags from a .filter file
		/// </summary>
		public static void LoadFilterFile(string filePath, out bool[] passedFilters)
		{
			using (FileStream inputStream = new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.Read))
			using (BinaryReader reader = new BinaryReader(inputStream))
			{
				int int1 = reader.ReadInt32();
				if (int1 == 0 && inputStream.Length == 4)
				{
					// This is a valid version-1 filter file with 0 clusters
					passedFilters = new bool[0];
					return;
				}

				int int2 = reader.ReadInt32();
				int clusterCount;
				if (int1 == 0 && int2 == 2) // obsolete filter file format
				{
					clusterCount = reader.ReadInt32();
					ushort[] filterData = new ushort[clusterCount];
					try
					{
						for (int clusterIndex = 0; clusterIndex < clusterCount; clusterIndex++)
						{
							filterData[clusterIndex] = reader.ReadUInt16();
						}
					}
					catch (EndOfStreamException)
					{
					}
					passedFilters = new bool[clusterCount];

					for (int clusterIndex = 0; clusterIndex < clusterCount; clusterIndex++)
					{
						passedFilters[clusterIndex] = ((filterData[clusterIndex] & 1) != 0);
					}
				}
				else
				{
					if (int1 == 0 && int2 != 3)
					{
						throw new ApplicationException("Unrecognized filter file id tag.");
					}
					if (int1 != 0)
					{
						// reset reader to beginning, and then parse old format filter files
						reader.BaseStream.Seek(0, SeekOrigin.Begin);
					}
					clusterCount = reader.ReadInt32();
					passedFilters = new bool[clusterCount];
					byte[] buffer = reader.ReadBytes(clusterCount);
					for (int clusterIndex = 0; clusterIndex < clusterCount; clusterIndex++)
					{
						if (buffer[clusterIndex] != 0) passedFilters[clusterIndex] = true;
					}
				}
			}
		}

		public static void SaveDemultiplexFile(ushort[] sampleNumbers, string filePath)
		{
			SaveDemultiplexFile(sampleNumbers, filePath, 0, sampleNumbers.Length);
		}

		public static void SaveDemultiplexFile(ushort[] sampleNumbers, string filePath, int StartRange, int PointCount)
		{
			using (FileStream outputStream = new FileStream(filePath, FileMode.Create))
			using (BinaryWriter writer = new BinaryWriter(outputStream))
			{
				const int version = 1;
				writer.Write(version);
				writer.Write(PointCount);
				for (int Index = StartRange; Index < StartRange + PointCount; Index++)
				{
					writer.Write(sampleNumbers[Index]);
				}
			}
		}

		public static ushort[] LoadDemultiplexFile(string filePath)
		{
			using (FileStream inputStream = new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.Read))
			using (BinaryReader reader = new BinaryReader(inputStream))
			{
				int version = reader.ReadInt32();
				if (version > 1)
				{
					throw new ApplicationException(string.Format("Unknown version number {0} in file {1}", version, filePath));
				}
				int count = reader.ReadInt32();
				ushort[] results = new ushort[count];
				byte[] buffer = reader.ReadBytes(count * sizeof(ushort));
				Buffer.BlockCopy(buffer, 0, results, 0, buffer.Length);

				return results;
			}
		}

		/// <summary>
		///     This utility function deals with the filter file naming mess.  New filter files are named
		///     like this: "s_1_0240.filter"  But along the way, RTA used other naming conventions.
		/// </summary>
		public static string GetFilterFilePath(int laneNumber, int tileNumber, string laneFolder)
		{
			string filterPath = Path.Combine(laneFolder, string.Format("s_{0}_{1:D4}.filter", laneNumber, tileNumber));
			if (!File.Exists(filterPath))
				filterPath = Path.Combine(laneFolder, string.Format("s_{0}_{1}.filter", laneNumber, tileNumber));
			int readNumber = 1;
			if (!File.Exists(filterPath))
				filterPath = Path.Combine(laneFolder, string.Format("s_{0}_{1}_{2:D4}.filter", laneNumber, readNumber, tileNumber));
			if (!File.Exists(filterPath))
				filterPath = Path.Combine(laneFolder, string.Format("s_{0}_{1}_{2}.filter", laneNumber, readNumber, tileNumber));
			string parentFolder = Path.GetDirectoryName(laneFolder);
			if (!File.Exists(filterPath))
				filterPath = Path.Combine(parentFolder, string.Format("s_{0}_{1:D4}.filter", laneNumber, tileNumber));
			return filterPath;
		}

#if MONO
    // TODO: Make these faster (native file reading/writing is a bit flakey under mono)
    
        public static void FastWrite(FileStream Stream, int[] Values)
        {
            byte[] buffer = new byte[Values.Length * sizeof(int)];
            Buffer.BlockCopy(Values, 0, buffer, 0, buffer.Length);
            Stream.Write(buffer, 0, buffer.Length);
        }

        public static void FastWrite(FileStream Stream, short[] Values)
        {
            byte[] buffer = new byte[Values.Length * sizeof(short)];
            Buffer.BlockCopy(Values, 0, buffer, 0, buffer.Length);
            Stream.Write(buffer, 0, buffer.Length);
        }

        public static void FastWrite(FileStream Stream, ushort[] Values)
        {
            byte[] buffer = new byte[Values.Length * sizeof(ushort)];
            Buffer.BlockCopy(Values, 0, buffer, 0, buffer.Length);
            Stream.Write(buffer, 0, buffer.Length);
        }

        public static void FastWrite(FileStream Stream, float[] Values)
        {
            byte[] buffer = new byte[Values.Length * sizeof(float)];
            Buffer.BlockCopy(Values, 0, buffer, 0, buffer.Length);
            Stream.Write(buffer, 0, buffer.Length);
        }
        
        public static float[] FastReadFloat(FileStream Stream, int Count)
        {
            float[] Results = new float[Count];
            byte[] buffer = new byte[Count * sizeof(float)];
            Stream.Read(buffer, 0, buffer.Length);
            Buffer.BlockCopy(buffer, 0, Results, 0, buffer.Length);
        
            return Results;
        }

        public static int[] FastReadInt(FileStream Stream, int Count)
        {
            int[] Results = new int[Count];
            byte[] buffer = new byte[Count * sizeof(int)];
            Stream.Read(buffer, 0, buffer.Length);
            Buffer.BlockCopy(buffer, 0, Results, 0, buffer.Length);
        
            return Results;
        }

        public static short[] FastReadShort(FileStream Stream, int Count)
        {
            short[] Results = new short[Count];
            byte[] buffer = new byte[Count * sizeof(short)];
            Stream.Read(buffer, 0, buffer.Length);
            Buffer.BlockCopy(buffer, 0, Results, 0, buffer.Length);
        
            return Results;
        }
#else
		public static unsafe void FastWrite(FileStream stream, int[] values)
		{
			fixed (int* data = values)
			{
				Win32Write(stream.SafeFileHandle.DangerousGetHandle(), (byte*)data, values.Length * sizeof(int));
			}
		}

		public static unsafe void FastWrite(FileStream stream, short[] values)
		{
			fixed (short* data = values)
			{
				Win32Write(stream.SafeFileHandle.DangerousGetHandle(), (byte*)data, values.Length * sizeof(short));
			}
		}

		public static unsafe void FastWrite(FileStream stream, ushort[] values)
		{
			fixed (ushort* data = values)
			{
				Win32Write(stream.SafeFileHandle.DangerousGetHandle(), (byte*)data, values.Length * sizeof(ushort));
			}
		}

		public static unsafe void FastWrite(FileStream stream, float[] values)
		{
			fixed (float* data = values)
			{
				Win32Write(stream.SafeFileHandle.DangerousGetHandle(), (byte*)data, values.Length * sizeof(float));
			}
		}

		public static unsafe float[] FastReadFloat(FileStream stream, int count)
		{
			float[] result = new float[count];
			fixed (float* data = result)
			{
				Win32Read(stream.SafeFileHandle.DangerousGetHandle(), (byte*)data, count * sizeof(float));
			}
			return result;
		}

		public static unsafe int[] FastReadInt(FileStream stream, int count)
		{
			int[] result = new int[count];
			fixed (int* data = result)
			{
				Win32Read(stream.SafeFileHandle.DangerousGetHandle(), (byte*)data, count * sizeof(int));
			}
			return result;
		}

		public static unsafe ushort[] FastReadUshort(FileStream stream, int count)
		{
			ushort[] result = new ushort[count];

			if (CrossPlatform.IsThisMono())
			{
				byte[] buffer = new byte[count * sizeof(ushort)];
				stream.Read(buffer, 0, buffer.Length);
				Buffer.BlockCopy(buffer, 0, result, 0, buffer.Length);
			}
			else
			{
				fixed (ushort* data = result)
				{
					Win32Read(stream.SafeFileHandle.DangerousGetHandle(), (byte*)data, count * sizeof(ushort));
				}
			}
			return result;
		}

		public static unsafe short[] FastReadShort(FileStream stream, int count)
		{
			short[] result = new short[count];
			fixed (short* data = result)
			{
				Win32Read(stream.SafeFileHandle.DangerousGetHandle(), (byte*)data, count * sizeof(ushort));
			}
			return result;
		}

		/// <summary>
		///     Unmanaged read.  Faster (marginally) than doing BinaryReader.ReadBytes() followed by
		///     System.Buffer.BlockCopy(), because we don't need to allocate two different arrays and copy
		///     from one to the other.
		/// </summary>
		private static unsafe void Win32Read(IntPtr handle, byte* buffPtr, long length)
		{
			UInt32 totalBytesRead = 0;
			UInt32 chunkSize = FileChunkSize;
			long remainingBytesToRead = length;
			bool result = true;
			while (remainingBytesToRead > 0)
			{
				result = true;
				UInt32 bytesRead;
				result = Functions.ReadFile(handle, (IntPtr)buffPtr,
											(remainingBytesToRead > chunkSize)
												? chunkSize
												: (UInt32)remainingBytesToRead,
											out bytesRead, IntPtr.Zero);

				if (!result)
				{
					chunkSize -= Megabyte;
					if (chunkSize < Megabyte) break;
				}
				else
				{
					buffPtr += bytesRead;
					totalBytesRead += bytesRead;
					remainingBytesToRead -= bytesRead;
				}
			}
			if (!result)
			{
				throw new IOException(string.Format("Failed the Win32 ReadFile() function, error code = {0}",
													Functions.GetLastError()));
			}
			if (totalBytesRead != length)
			{
				throw new IOException(string.Format(
					"Win32 WriteFile() function failed to read the desired number of bytes ({0} requested, {1} written)",
					length, totalBytesRead));
			}
		}


		private const int Megabyte = 1024 * 1024;
		private const int FileChunkSize = 50 * Megabyte; // 50 Mb

		/// <summary>
		///     Unmanaged write, for speed.
		/// </summary>
		private static unsafe void Win32Write(IntPtr handle, byte* buffPtr, long length)
		{
			UInt32 totalBytesWritten = 0;

			// Write file header from memory to file in one fast write
			bool result = true;
			long remainBytesToWrite = length;
			UInt32 chunkSize = FileChunkSize;

			while (remainBytesToWrite > 0)
			{
				result = true;
				UInt32 bytesWritten;
				result = Functions.WriteFile(handle, (IntPtr)buffPtr,
											 (remainBytesToWrite > chunkSize) ? chunkSize : (UInt32)remainBytesToWrite,
											 out bytesWritten, IntPtr.Zero);
				// We failed to write using the specified chunk size, decrement by a megabyte and try again.
				if (!result)
				{
					chunkSize -= Megabyte;

					// If the chunk size is less than a megabyte then we need to throw an exception.
					if (chunkSize < Megabyte)
						break;
				}
				else
				{
					// Increment the buffer
					buffPtr += bytesWritten;
					totalBytesWritten += bytesWritten;
					remainBytesToWrite -= bytesWritten;
				}
			}

			if (!result)
			{
				throw new IOException(string.Format(
					"Failed the Win32 WriteFile() function, error code = {0}",
					Functions.GetLastError()));
			}

			if (totalBytesWritten != length)
			{
				throw new IOException(string.Format(
					"Win32 WriteFile() function failed to write the desired number of bytes ({0} requested, {1} written)",
					length, totalBytesWritten));
			}
		}


#endif
	}
}