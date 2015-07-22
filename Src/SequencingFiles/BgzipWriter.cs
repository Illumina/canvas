using System;
using System.Text;
using System.IO;
using System.Collections.Generic;

namespace SequencingFiles
{
	/// <summary>
	///     Writes to bgzip format
	///     Note there is no reader at this time because any valid bgzip file can be read with a gzip reader,
	///     however we may need some read functionality if we bring vcf indexing into C#.
	/// </summary>
	public class BgzipWriter : BgzfWriterCommon, BgzipOrStreamWriter.IWriter
	{
		#region member variables

		private const byte Newline = 10;

		private const int InitBufSize = 4 * 1024;
		private readonly ASCIIEncoding _encoding = new ASCIIEncoding();
		private byte[] _outputBuffer = new byte[InitBufSize];

		#endregion

		public BgzipWriter(string filename, GzipCompressionLevel compressionLevel = GzipCompressionLevel.DefaultCompression5)
			: base((int)compressionLevel)
		{
			base.Open(filename);
		}

		public BgzipWriter(Stream outStream, GzipCompressionLevel compressionLevel = GzipCompressionLevel.DefaultCompression5)
			: base((int)compressionLevel)
		{
			base.Open(outStream);
		}

		public void Write(byte[] buffer)
		{
			this.Write(buffer, (uint)buffer.Length);
		}

		/// <summary>
		///     Writes a line to the file
		/// </summary>
		public void WriteLine(string s = "")
		{
			int size = (_encoding.GetByteCount(s) + 1);
			if (size > _outputBuffer.Length)
			{
				_outputBuffer = new byte[size + 1024];
			}

			int outputSize = _encoding.GetBytes(s, 0, s.Length, _outputBuffer, 0);
			if ((outputSize + 1) != size)
			{
				throw new Exception(String.Format("ERROR: failed to byte-encode string: \"{0}\"", s));
			}
			_outputBuffer[outputSize++] = Newline;
			Write(_outputBuffer, (uint)outputSize);
		}

		/// <summary>
		///     Writes to the file without a newline
		/// </summary>
		public void Write(string s)
		{
			int size = _encoding.GetByteCount(s);
			if (size > _outputBuffer.Length)
			{
				_outputBuffer = new byte[size + 1024];
			}

			int outputSize = _encoding.GetBytes(s, 0, s.Length, _outputBuffer, 0);
			if (outputSize != size)
			{
				throw new Exception(String.Format("ERROR: failed to byte-encode string: \"{0}\"", s));
			}
			Write(_outputBuffer, (uint)outputSize);
		}

		public void WriteLine(string format, params Object[] arg)
		{
			WriteLine(string.Format(format, arg));
		}

		public void Write(string format, params Object[] arg)
		{
			Write(string.Format(format, arg));
		}

		/// <summary>
		///     Bgzip a text file. Also converts any \r\n newlines to \n
		/// </summary>
		public static void CompressFile(string inputFilename, string outputFilename, GzipCompressionLevel compressionLevel = GzipCompressionLevel.DefaultCompression5)
		{
			using (FileStream inStream = new FileStream(inputFilename, FileMode.Open, FileAccess.Read, FileShare.Read))
			using (StreamReader reader = new StreamReader(inStream))
			{
				CompressFile(reader, new FileStream(outputFilename, FileMode.Create), compressionLevel);
			}
		}

		/// <summary>
		///     Bgzip a text file. Also converts any \r\n newlines to \n
		/// </summary>
		public static void CompressFile(StreamReader inputStream, Stream outputStream, GzipCompressionLevel compressionLevel = GzipCompressionLevel.DefaultCompression5)
		{
			// our output file
			using (BgzipWriter writer = new BgzipWriter(outputStream, compressionLevel))
			{
				string fileLine;
				while ((fileLine = inputStream.ReadLine()) != null)
				{
					writer.WriteLine(fileLine);
				}
			}
		}

		/// <summary>
		///     Given a list of input bgzipped files, output a single bgzipped file which is a concatentation
		///     of the input files. The zero size blocks that are used as EOF markers will be removed except for a final 
		///     zero size block at the end of the output file.
		///     Note: this was adapted from bgzf_cat.c written by Chris Saunders
		/// </summary>
		public static void BgzfCat(List<string> inputFilenames, string outputFilename, GzipCompressionLevel compressionLevel = GzipCompressionLevel.DefaultCompression5)
		{
			// some useful constants:
			int gZipId1 = 31;
			int gZipId2 = 139;
			int emptyBlockSize = 28;

			// our buffers
			byte[] blockBuffer = new byte[MaxBlockSize];
			byte[] emptyBlockBuffer = new byte[emptyBlockSize];

			// our output file
			using (BgzipWriter writer = new BgzipWriter(outputFilename, compressionLevel))
			{
				foreach (string inputFilename in inputFilenames)
				{
					using (BgzipReader inputFile = new BgzipReader(inputFilename))
					{
						// if we uncompressed some data in the process of opening this file we need to compress and write it back to the output
						if (inputFile.BlockOffset < inputFile.BlockLength)
						{
							//sanity check, we just opened this file so the offset should be 0
							if (inputFile.BlockOffset != 0)
								throw new ApplicationException("Opened bgzipped file and not at the start of a new block!");
							writer.Write(inputFile.LineBuffer, (uint)(inputFile.BlockLength - inputFile.BlockOffset));
							writer.FlushBlock();
						}
						bool firstRead = true;
						int bytesRead;
						while ((bytesRead = inputFile.BgzfFileStream.Read(blockBuffer, 0, MaxBlockSize)) > 0)
						{
							if (bytesRead < emptyBlockSize)
							{
								// if this is the first read it must be at least the size of an empty block
								if (firstRead)
									throw new ApplicationException("BgzfCat Error: truncated file?: " + inputFilename);

								// this is the remainder of the final empty block for this input file

								// let's write out the beginning of the data from the emptyBlockBuffer which we know is not part of this final empty block (i.e. it is the end of the block from the previous write)
								writer.BgzfFileStream.Write(emptyBlockBuffer, 0, bytesRead);

								// save the final empty block so we can perform a sanity check
								Buffer.BlockCopy(emptyBlockBuffer, bytesRead, emptyBlockBuffer, 0, emptyBlockSize - bytesRead);
								Buffer.BlockCopy(blockBuffer, 0, emptyBlockBuffer, emptyBlockSize - bytesRead, bytesRead);
							}
							else
							{
								if (!firstRead)
								{
									// let's write the data that we saved from the previous read (i.e. it wasn't the empty block EOF)
									writer.BgzfFileStream.Write(emptyBlockBuffer, 0, emptyBlockSize);
								}
								// save some data at the end in case it is from the empty block EOF. we don't want to write out the empty block EOF.
								// the empty block EOF is automatically written when the BgzipWriter is disposed
								Buffer.BlockCopy(blockBuffer, bytesRead - emptyBlockSize, emptyBlockBuffer, 0, emptyBlockSize);
								//write out the data that we aren't saving
								writer.BgzfFileStream.Write(blockBuffer, 0, bytesRead - emptyBlockSize);
							}
							firstRead = false;
						}

						// if the one and only block was the empty block we need to copy it to emptyBlockBuffer to perform sanity checking
						if (firstRead)
							Buffer.BlockCopy(inputFile.CompressedBlock, 0, emptyBlockBuffer, 0, emptyBlockSize);

						// sanity check for the final gzip block 
						int blockSize = emptyBlockBuffer[emptyBlockSize - 4];
						if (emptyBlockBuffer[0] != gZipId1 || emptyBlockBuffer[1] != gZipId2 || blockSize != 0)
						{
							throw new ApplicationException("BgzfCat Error: unexpected final block structure in file " + inputFilename);
						}
					}
				}
			}
		}
	}


	/// <summary>
	/// Simple Wrapper that can write to either a StreamWriter (plain-text) or a BgzipWriter (gzipped output)
	/// </summary>
	public class BgzipOrStreamWriter : ClosableDisposable
	{
		internal interface IWriter : IClosableDisposable
		{
			void Write(string s);
			void Write(string s, params object[] args);
			void WriteLine(string s);
			void WriteLine(string s, params object[] args);
		}

		private class InternalStreamWriter : StreamWriter, IWriter
		{
			public InternalStreamWriter(string filename) : base(filename) { }
		}
		IWriter writer;

		/// <summary>
		/// Create a new Writer writing to filename. If filename ends in .gz output will be bgzip compressed, otherwise plain-text
		/// </summary>
		public BgzipOrStreamWriter(string filename, GzipCompressionLevel compression = GzipCompressionLevel.DefaultCompression5)
		{
			if (filename.ToLower().EndsWith(".gz"))
			{
				writer = new BgzipWriter(filename, compression);
			}
			else
			{
				var streamWriter = new InternalStreamWriter(filename);
				streamWriter.NewLine = "\n";
				writer = streamWriter;
			}
		}

		public void Write(string s)
		{
			writer.Write(s);
		}

		public void Write(string s, params object[] args)
		{
			writer.Write(s, args);
		}

		public void WriteLine(string s = "")
		{
			writer.WriteLine(s);
		}

		public void WriteLine(string s, params object[] args)
		{
			writer.WriteLine(s, args);
		}

		public override void Close()
		{
			writer.Close();
		}
	}
}