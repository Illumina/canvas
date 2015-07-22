using System;
using System.Collections.Generic;
using System.IO;
using SequencingFiles.Compression;

namespace SequencingFiles
{
	public class GzipReader : GzipCommon
	{
		#region member variables
		private const uint BufferSize = 131072;
		private readonly char[] _stringBuffer;
		private int _bufferByteCount;
		#endregion

		/// <summary>
		///     constructor
		/// </summary>
		public GzipReader(string filename)
		{
			LineBuffer = new byte[BufferSize];
			_stringBuffer = new char[BufferSize];
			_bufferByteCount = 0;
			Open(filename);
		}

		/// <summary>
		/// uncompress foo.gz to foo
		/// this will probably fail if there are no newlines in the file and the entire file cannot fit into memory
		/// TODO: implement GzipReader.ReadBytes so that we don't need to use ReadLine/WriteLine
		/// </summary>
		public static void UncompressFile(string sourcePath, string targetPath)
		{
			using (StreamWriter writer = new StreamWriter(targetPath))
			using (GzipReader reader = new GzipReader(sourcePath))
			{
				writer.NewLine = "\n";
				while (true)
				{
					string fileLine = reader.ReadLine();
					if (fileLine == null) break;
					writer.WriteLine(fileLine);
				}
			}
		}

		/// <summary>
		///     Closes the FASTQ file
		/// </summary>
		public override void Close()
		{
			if (IsOpen)
			{
				IsOpen = false;
				SafeNativeMethods.gzclose(FileStreamPointer);
			}
		}

		/// <summary>
		///     Fills the buffer
		/// </summary>
		/// <param name="bufferOffset"></param>
		public void FillBuffer(int bufferOffset)
		{
			_bufferByteCount = SafeNativeMethods.gzreadOffset(FileStreamPointer, LineBuffer, bufferOffset,
															 BufferSize - (uint)bufferOffset);
			if (_bufferByteCount < 0)
			{
				throw new ApplicationException(string.Format("ERROR: Unable to read data from {0}, _bufferByteCount={1}, check zlib.h or zutil.c for error code meaning.", FilePath,
					_bufferByteCount));
			}
		}

		/// <summary>
		///     returns a null-terminated string from a binary reader
		/// </summary>
		private static string GetNullTerminatedString(BinaryReader reader)
		{
			List<char> nameBytes = new List<char>();

			while (true)
			{
				byte b = reader.ReadByte();
				if (b == 0) break;
				nameBytes.Add((char)b);
			}

			return new string(nameBytes.ToArray());
		}

		/// <summary>
		///     returns the uncompressed file size of a compressed file.
		///     N.B. Only works for bgzipped files at the moment
		///     returns -1 if the file size cannot be determined.
		/// </summary>
		public static long GetUncompressedSize(string filePath)
		{
			long uncompressedSize = 0;

			using (FileStream fileStream = new FileStream(filePath, FileMode.Open))
			using (BinaryReader reader = new BinaryReader(fileStream))
			{
				try
				{
					while (true)
					{
						uncompressedSize += GetUncompressedSizeFromGzipHeader(reader);
					}
				}
				catch (EndOfStreamException)
				{
					// this how EOFs are detected in binary files... lame.
				}
				catch (Exception)
				{
					uncompressedSize = -1;
				}
			}

			return uncompressedSize;
		}

		/// <summary>
		/// parses the header of a gzip/bgzip file and determines the
		/// uncompressed file size in that compression block.
		/// </summary>
		private static long GetUncompressedSizeFromGzipHeader(BinaryReader reader)
		{
			// make sure we get the magic bytes 31 & 139 (=35615)
			ushort id = reader.ReadUInt16();
			if (id != 35615) throw new ApplicationException("Found a bad header.");

			reader.ReadByte();
			byte flags = reader.ReadByte();
			reader.ReadUInt32();
			reader.ReadUInt16();

			int numHeaderBytes = 10;

			// evaluate our fields
			//bool hasText        = (flags & 1)  != 0;
			bool hasCrc = (flags & 2) != 0;
			bool hasExtraFields = (flags & 4) != 0;
			bool hasName = (flags & 8) != 0;
			bool hasComment = (flags & 16) != 0;

			// currently we only support bgzip files - so extra fields are required
			if (!hasExtraFields) throw new ApplicationException("Missing extra fields.");

			// load the extra fields
			reader.ReadUInt16();
			ushort subfieldId = reader.ReadUInt16();

			// sanity check: make sure this file uses BGZF blocks
			if (subfieldId != 17218) throw new ApplicationException("Could not find the BGZF subfield IDs.");

			reader.ReadUInt16();
			ushort blockSize = reader.ReadUInt16();

			numHeaderBytes += 8;

			// load the name
			if (hasName)
			{
				string name = GetNullTerminatedString(reader);
				numHeaderBytes += name.Length + 1;
			}

			// load the comment
			if (hasComment)
			{
				string comment = GetNullTerminatedString(reader);
				numHeaderBytes += comment.Length + 1;
			}

			// load the CRC
			if (hasCrc)
			{
				reader.ReadUInt16();
				numHeaderBytes += 2;
			}

			// skip some bytes
			numHeaderBytes += 8;
			long compressedSize = blockSize - numHeaderBytes + 1;
			reader.BaseStream.Seek(compressedSize, SeekOrigin.Current);

			// read the crc and uncompressed file size
			// BGZF blocks are typically 64K, so we're fine with the 32-bit variable
			reader.ReadUInt32();
			uint uncompressedSize = reader.ReadUInt32();

			return uncompressedSize;
		}

		/// <summary>
		///     Gets the next string in the file.
		/// </summary>
		/// <returns>Returns the next string or null.</returns>
		public string ReadLine()
		{
			string s = "";

			if (_bufferByteCount < 0)
			{
				//throw exception if this is <0
				throw new ApplicationException(string.Format("ERROR: Unable to read data from {0}, _bufferByteCount={1}, check zlib.h or zutil.c for error code meaning.", FilePath,
					_bufferByteCount));
			}
			// skip if the file is not currently open or if we don't have any data in the buffer
			if (!IsOpen || (_bufferByteCount <= 0)) return null;


			int crOffset = -1;
			while (true)
			{

				crOffset = Array.IndexOf(LineBuffer, LineFeedChar, CurrentOffset, _bufferByteCount - CurrentOffset);

				if (crOffset != -1)
				{
					s = s + GetString(crOffset - CurrentOffset);
					CurrentOffset = crOffset + 1;
					break;
				}
				else
				{
					int remainingLen = _bufferByteCount - CurrentOffset;   // remain of this 
					s = s + GetString(remainingLen);

					FillBuffer(0);

					if (_bufferByteCount > 0)
					{
						BufferOffset += CurrentOffset + remainingLen;
						CurrentOffset = 0;
					}
					else if (_bufferByteCount == 0)
					{
						return string.IsNullOrEmpty(s) ? null : s;
					}
					else
					{
						//throw exception if _bufferByteCount <0
						throw new ApplicationException(string.Format("ERROR: Unable to read data from {0}, _bufferByteCount={1}, check zlib.h or zutil.c for error code meaning.", FilePath,
							_bufferByteCount));
					}
				}
			}

			return s;
		}

		/// <summary>
		///     Converts a byte array into a string
		/// </summary>
		/// <param name="len"></param>
		/// <returns></returns>
		private string GetString(int len)
		{
			for (int charIndex = 0; charIndex < len; charIndex++)
			{
				_stringBuffer[charIndex] = (char)LineBuffer[CurrentOffset + charIndex];
			}
			return new string(_stringBuffer, 0, len);
		}

		/// <summary>
		///     Opens the file
		/// </summary>
		public void Open(string filename)
		{
			Open(filename, "rb");
			FillBuffer(0);
		}

		/// <summary>
		///     Moves the file stream pointer to the beginning of the file
		/// </summary>
		public void Rewind()
		{
			SafeNativeMethods.gzrewind(FileStreamPointer);
			CurrentOffset = 0;
			BufferOffset = 0;
			FillBuffer(0);
		}
	}
}