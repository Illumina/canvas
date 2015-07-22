using SequencingFiles.Compression;
using System.IO;

namespace SequencingFiles
{
	/// <summary>
	///     Defines the possible GZIP compression levels
	/// </summary>
	public enum GzipCompressionLevel
	{
		BestSpeed1 = 1,
		CompressionLevel2 = 2,
		CompressionLevel3 = 3,
		CompressionLevel4 = 4,
		DefaultCompression5 = 5,
		CompressionLevel6 = 6,
		CompressionLevel7 = 7,
		CompressionLevel8 = 8,
		BestCompression9 = 9
	}

	public class GzipWriter : GzipCommon
	{
		#region member variables
		private const uint GwBufferSize = 65536;
		#endregion

		/// <summary>
		///     constructor
		/// </summary>
		public GzipWriter(string filename, GzipCompressionLevel compLevel = GzipCompressionLevel.DefaultCompression5, bool append = false)
		{
			LineBuffer = new byte[GwBufferSize];
			Open(filename, compLevel, append);
		}

		public static void Uncompress(string tempPath, string outputPath)
		{
			if (!File.Exists(tempPath))
			{
				return;
			}
			using (GzipReader reader = new GzipReader(tempPath))
			using (StreamWriter writer = new StreamWriter(outputPath))
			{
				writer.NewLine = "\n";
				while (true)
				{
					string FileLine = reader.ReadLine();
					if (FileLine == null) break;
					writer.WriteLine(FileLine);
				}
			}
		}

		/// <summary>
		///     Adds the string to the output buffer
		/// </summary>
		/// <param name="charArray"></param>
		/// <param name="offset"></param>
		/// <param name="len"></param>
		private void AddString(ref char[] charArray, int offset, int len)
		{
			for (int charIndex = 0; charIndex < len; charIndex++)
			{
				LineBuffer[CurrentOffset++] = (byte)charArray[offset + charIndex];
			}
		}

		/// <summary>
		///     Closes the file
		/// </summary>
		public override void Close()
		{
			if (IsOpen)
			{
				IsOpen = false;
				FlushBuffer();
				SafeNativeMethods.gzclose(FileStreamPointer);
			}
		}

		/// <summary>
		///     Flushes the output buffer to disk
		/// </summary>
		private void FlushBuffer()
		{
			SafeNativeMethods.gzwriteOffset(FileStreamPointer, LineBuffer, 0, (uint)CurrentOffset);
			CurrentOffset = 0;
		}

		/// <summary>
		///     Opens the file
		/// </summary>
		/// <param name="filename">the file to open</param>
		/// <param name="compLevel">the compression level</param>
		/// <param name="appendFlag"></param>
		public void Open(string filename, GzipCompressionLevel compLevel, bool appendFlag)
		{
			string fileMode = string.Format("wb{0}", (int)compLevel);
			if (appendFlag) fileMode = string.Format("ab{0}", (int)compLevel);
			base.Open(filename, fileMode);
		}

		/// <summary>
		///     Writes a line to the file
		/// </summary>
		public void WriteLine(string s)
		{
			// skip if the file is not currently open
			if (!IsOpen) return;

			int remainingBytes = s.Length;
			int availableBufferBytes = (int)GwBufferSize - CurrentOffset;
			int stringOffset = 0;

			char[] charArray = s.ToCharArray();

			while (remainingBytes > availableBufferBytes)
			{
				// fill the buffer
				AddString(ref charArray, stringOffset, availableBufferBytes);

				// flush the buffer
				FlushBuffer();

				// update
				stringOffset += availableBufferBytes;
				remainingBytes -= availableBufferBytes;

				availableBufferBytes = (int)GwBufferSize - CurrentOffset;
			}

			// add the remaining bytes to the buffer
			if (remainingBytes > 0)
			{
				AddString(ref charArray, stringOffset, remainingBytes);
				availableBufferBytes = (int)GwBufferSize - CurrentOffset;
			}

			// add a carriage return
			if (availableBufferBytes == 0) FlushBuffer();
			LineBuffer[CurrentOffset++] = LineFeedChar;
		}
	}
}