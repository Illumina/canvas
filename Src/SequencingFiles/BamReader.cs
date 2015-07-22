using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using SequencingFiles.Compression;

namespace SequencingFiles
{
	public interface IBamReader : IDisposable
	{
		bool GetNextAlignment(ref BamAlignment alignment, bool skipAdditionalParsing);
		int GetReferenceIndex(string referenceName);
		bool Jump(int refID, int position);
		List<string> GetReferenceNames();
	}

	public class BamReader : BgzfCommon, IBamReader
	{
		#region member variables

		private long _alignmentsOffset;
		private char[] _baseLookupTable;
		private byte[] _byteBuffer;
		private bool _hasIndex;
		private string _header;

		private BamIndex _index;
		private BinaryReader _reader;
		private List<GenomeMetadata.SequenceMetadata> _referenceIndex;
		private Dictionary<string, int> _referenceNameToIndex;
		private char[] _sequenceBuffer;

		public string BamPath;

		#endregion

		// constructor
		public BamReader()
		{
			Initialize();
		}

		// constructor
		public BamReader(string filename)
		{
			Initialize();
			Open(filename);
		}

		// close BAM file
		public override void Close()
		{
			if (!IsOpen) return;
			IsOpen = false;

			_header = "";

			_reader.Close();
			BgzfFileStream.Close();
		}

		// returns the SAM header
		public string GetHeader()
		{
			return _header;
		}

		/// <summary>
		/// access UR fields from header text. ISAAC provides this, but bwa does not
		/// one list entry per sequence
		/// list entry is null if UR field does not exist in header
		/// </summary>
		public List<string> GetSequenceURIs()
		{
			List<string> uris = new List<string>();
			string[] lines = _header.Split('\n');
			foreach (string line in lines)
			{
				if (!line.StartsWith("@SQ")) continue;
				string[] fields = line.Split('\t');
				string field = fields.FirstOrDefault(s => s.StartsWith("UR:"));
				if (field != null)
					field = field.Substring(3).Trim();
				uris.Add(field);
			}
			return uris;
		}

		/// <summary>
		/// Return SM tags, one for each read group (@RG). 
		/// List entry is null if SM is missing for corresponding read group.
		/// </summary>
		public List<string> GetReadGroupSamples()
		{
			List<string> samples = new List<string>();
			string[] lines = _header.Split('\n');
			foreach (string line in lines)
			{
				if (!line.StartsWith("@RG")) continue;
				string[] fields = line.Split('\t');
				string field = fields.FirstOrDefault(s => s.StartsWith("SM:"));
				if (field != null)
					field = field.Substring(3).Trim();
				samples.Add(field);
			}
			return samples;
		}

		/// <summary>
		/// return first non-null non-empty uri from the header
		/// return null if we can't find one (e.g. bwa does not provide UR field)
		/// </summary>
		/// <returns></returns>
		public string GetReferenceURI()
		{
			return GetSequenceURIs().FirstOrDefault(s => !string.IsNullOrEmpty(s));
		}

		/// <summary>
		///     Returns the remainder of the SAM header (reference sequence names and lengths):
		/// </summary>
		public string GetSequenceHeaderString()
		{
			StringBuilder builder = new StringBuilder();
			foreach (GenomeMetadata.SequenceMetadata refSeq in _referenceIndex)
			{
				builder.AppendFormat("@SQ\tSN:{0}\tLN:{1}\n", refSeq.Name, refSeq.Length);
			}
			return builder.ToString();
		}

		// retrieves next available alignment
		public bool GetNextAlignment(ref BamAlignment alignment, bool skipAdditionalParsing)
		{
			// check that our file is open
			if (!IsOpen) return false;

			// retrieve the alignment data length
			if (Read(ref _byteBuffer, 4) != 4) return false;
			uint alignmentDataLen = BitConverter.ToUInt32(_byteBuffer, 0);
			if (alignmentDataLen == 0) return false;

			// retrieve the alignment data
			if (Read(ref _byteBuffer, alignmentDataLen) != alignmentDataLen) return false;

			// retrieve the core alignment data
			uint compositeData1 = BitConverter.ToUInt32(_byteBuffer, 8);
			uint flagAndNumCigarOps = BitConverter.ToUInt32(_byteBuffer, 12);
			uint numBases = BitConverter.ToUInt32(_byteBuffer, 16);
			if (numBases > _sequenceBuffer.Length)
			{
				// For very long reads, re-allocate this buffer to twice the data length
				_sequenceBuffer = new char[numBases * 2];
			}
			uint readNameLen = compositeData1 & 0xff;
			uint numCigarOps = flagAndNumCigarOps & 0xffff;

			alignment.RefID = BitConverter.ToInt32(_byteBuffer, 0);
			alignment.Position = BitConverter.ToInt32(_byteBuffer, 4);
			alignment.Bin = (compositeData1 >> 16);
			alignment.MapQuality = ((compositeData1 >> 8) & 0xff);
			alignment.AlignmentFlag = flagAndNumCigarOps >> 16;
			alignment.MateRefID = BitConverter.ToInt32(_byteBuffer, 20);
			alignment.MatePosition = BitConverter.ToInt32(_byteBuffer, 24);
			alignment.FragmentLength = BitConverter.ToInt32(_byteBuffer, 28);

			// retrieve the read name
			int offset = (int)BamConstants.CoreAlignmentDataLen;
			alignment.Name = Encoding.ASCII.GetString(_byteBuffer, offset, (int)(readNameLen - 1));
			offset += (int)readNameLen;

			// retrieve the CIGAR operations
			alignment.CigarData.Clear();
			for (uint i = 0; i < numCigarOps; ++i, offset += 4)
			{
				uint cigarData = BitConverter.ToUInt32(_byteBuffer, offset);
				alignment.CigarData.Add(new CigarOp(BamConstants.CigarTypes[cigarData & BamConstants.CigarMask],
													cigarData >> BamConstants.CigarShift));
			}

			// here we provide a mechanism for skipping the processing of
			// bases, base qualities, and tags
			if (!skipAdditionalParsing)
			{
				// retrieve the bases
				byte shift = 4;
				for (int i = 0; i < numBases; ++i, shift ^= 4)
				{
					_sequenceBuffer[i] = _baseLookupTable[(_byteBuffer[offset] >> shift) & 15];
					if (shift == 0) offset++;
				}

				if (shift == 0) offset++;

				alignment.Bases = new string(_sequenceBuffer, 0, (int)numBases);

				// retrieve the qualities
				if ((alignment.Qualities == null) || (alignment.Qualities.Length != numBases))
				{
					alignment.Qualities = new byte[numBases];
				}

				Buffer.BlockCopy(_byteBuffer, offset, alignment.Qualities, 0, (int)numBases);
				offset += (int)numBases;

				// retrieve the tags
				int numTagBytes = (int)alignmentDataLen - offset;
				alignment.TagData = new byte[numTagBytes];
				Array.Copy(_byteBuffer, offset, alignment.TagData, 0, numTagBytes);
			}

			return true;
		}

		/// <summary>
		///     given a reference name, this method returns the reference index
		/// </summary>
		public int GetReferenceIndex(string referenceName)
		{
			int ret;
			if (!_referenceNameToIndex.TryGetValue(referenceName, out ret)) ret = -1;
			return ret;
		}

		// returns the SAM header
		public List<GenomeMetadata.SequenceMetadata> GetReferences()
		{
			return _referenceIndex;
		}

		/// <summary>
		///     returns the reference sequence name found the given reference ID
		/// </summary>
		public string GetReferenceNameByID(int refID)
		{
			return _referenceIndex[refID].Name;
		}

		/// <summary>
		///     returns the reference sequence names found in this BAM file
		/// </summary>
		public List<string> GetReferenceNames()
		{
			return _referenceIndex.Select(sm => sm.Name).ToList();
		}

		// common initialization routines used by the constructors
		private void Initialize()
		{
			_byteBuffer = new byte[MaxBlockSize];
			_sequenceBuffer = new char[BamConstants.MaxReadLength];
			_baseLookupTable = "=ACMGRSVTWYHKDBN".ToCharArray();

			_referenceNameToIndex = new Dictionary<string, int>();
			_referenceIndex = new List<GenomeMetadata.SequenceMetadata>();
			_alignmentsOffset = 0;

			_index = new BamIndex();
		}

		/// <summary>
		///     returns true if the alignment overlaps with the specified interval
		/// </summary>
		private bool IsOverlap(int begin, int end, BamAlignment alignment)
		{
			int alignmentBegin = alignment.Position;
			int alignmentEnd = alignment.GetEndPosition();
			return (alignmentEnd >= begin) && (alignmentBegin < end);
		}

		/// <summary>
		///     jumps to the specified position in the BAM file
		/// </summary>
		/// <returns>true if we were successfully able to jump to the requested position</returns>
		public bool Jump(string referenceName, int position)
		{
			int referenceIndex;
			if (!_referenceNameToIndex.TryGetValue(referenceName, out referenceIndex))
			{
				throw new ApplicationException(
					string.Format("Unable to find the reference sequence ({0}) in the BAM file.", referenceName));
			}

			return Jump(referenceIndex, position);
		}

		/// <summary>
		/// Jump to unaligned reads (with no associated chromosome) at end of bam file.
		/// </summary>
		public bool JumpToUnaligned()
		{
			// sanity check: make sure we have unaligned reads
			if (_index.NumUnalignedWithoutCoordinates == 0) return false;

			// get the last indexed BAM offset
			ulong currentOffset = _index.GetLargestBamOffset();

			// reposition our BAM reader
			if (currentOffset != 0)
			{
				Seek(currentOffset);
			}
			else
			{
				Rewind();
				currentOffset = Tell();
			}

			// skip all of the alignments that are aligned
			BamAlignment alignment = new BamAlignment();

			while (true)
			{
				// look for the desired alignment
				if (GetNextAlignment(ref alignment, false))
				{
					if (alignment.RefID == -1) break;
					currentOffset = Tell();
				}
				else
				{
					//  end of file or error
					return false;
				}
			}

			// reset the file position (since we already read blew past the good alignment)
			Seek(currentOffset);

			return true;
		}

		/// <summary>
		///     jumps to the specified position in the BAM file
		/// </summary>
		/// <returns>true if we were successfully able to jump to the requested position</returns>
		public bool Jump(int refID, int position)
		{
			// sanity checks
			if (!_hasIndex) return false;
			if (refID > _referenceIndex.Count) return false;
			if (position > _referenceIndex[refID].Length) return false;

			// calculate the candidate index regions
			BamIterator bamIterator;
			bool foundOffset = _index.GetOffsets(refID, position, out bamIterator);

			if (!foundOffset || (bamIterator.Offsets == null) || (bamIterator.Offsets.Length == 0)) return false;

			int currentOffsetIndex = -1;
			int lastOffsetIndex = bamIterator.Offsets.Length - 1;
			BamAlignment alignment = new BamAlignment();

			while (true)
			{
				// jump to the next chunk
				if ((bamIterator.CurrentOffset == 0) ||
					(bamIterator.CurrentOffset >= bamIterator.Offsets[currentOffsetIndex].End))
				{
					// no more chunks
					if (currentOffsetIndex == lastOffsetIndex) return false;

					// sanity check
					if ((currentOffsetIndex >= 0) &&
						(bamIterator.CurrentOffset != bamIterator.Offsets[currentOffsetIndex].End))
					{
						throw new ApplicationException(
							string.Format(
								"Found a potential bug in the BAM index routines. CurrentOffset ({0}) != Offsets[currentOffsetIndex].End ({1}",
								bamIterator.CurrentOffset, bamIterator.Offsets[currentOffsetIndex].End));
					}

					// not adjacent chunks; then seek
					if ((currentOffsetIndex < 0) || (bamIterator.Offsets[currentOffsetIndex].End != bamIterator.Offsets[currentOffsetIndex + 1].Begin))
					{
						Seek(bamIterator.Offsets[currentOffsetIndex + 1].Begin);
						bamIterator.CurrentOffset = Tell();
					}

					currentOffsetIndex++;
				}

				// look for the desired alignment
				if (GetNextAlignment(ref alignment, false))
				{
					// no need to proceed
					if ((alignment.RefID != bamIterator.RefID) || (alignment.Position >= bamIterator.End))
					{
						return false;
					}
					if (IsOverlap(bamIterator.Begin, bamIterator.End, alignment))
					{
						// this is the read we're looking for
						break;
					}
					bamIterator.CurrentOffset = Tell();
				}
				else
				{
					//  end of file or error
					return false;
				}
			}

			// reset the file position (since we already read blew past the good alignment)
			Seek(bamIterator.CurrentOffset);

			return true;
		}

		// loads the header data
		private void LoadHeaderData()
		{
			// check to see if proper BAM header
			byte[] buffer = new byte[4];

			if (Read(ref buffer, 4) != 4)
			{
				throw new ApplicationException("ERROR: Could not read the BAM magic number.");
			}

			string magicNumberString = Encoding.ASCII.GetString(buffer, 0, 4);

			if (magicNumberString != BamConstants.MagicNumber)
			{
				throw new ApplicationException(
					string.Format("ERROR: Expected the BAM magic number to be {0}, but found {1}.",
								  BamConstants.MagicNumber, magicNumberString));
			}

			// get BAM header text length
			Read(ref buffer, 4);
			uint headerTextLength = BitConverter.ToUInt32(buffer, 0);

			// get BAM header text
			byte[] headerText = new byte[headerTextLength];
			Read(ref headerText, headerTextLength);
			_header = Encoding.ASCII.GetString(headerText, 0, (int)headerTextLength);
		}

		// loads the reference data
		private void LoadReferenceData()
		{
			// get number of reference sequences
			byte[] buffer = new byte[1024];
			Read(ref buffer, 4);
			uint numRefSeqs = BitConverter.ToUInt32(buffer, 0);
			if (numRefSeqs == 0) return;

			// populate the references list

			for (uint chromosomeIndex = 0; chromosomeIndex < numRefSeqs; chromosomeIndex++)
			{
				// retrieve the reference name
				Read(ref buffer, 4);
				uint nameLen = BitConverter.ToUInt32(buffer, 0);

				Read(ref buffer, nameLen);
				string name = Encoding.ASCII.GetString(buffer, 0, (int)(nameLen - 1));

				// retrieve the reference length
				Read(ref buffer, 4);
				uint refLen = BitConverter.ToUInt32(buffer, 0);

				_referenceIndex.Add(new GenomeMetadata.SequenceMetadata(name, refLen, (int)chromosomeIndex));
				_referenceNameToIndex[name] = (int)chromosomeIndex;
			}
		}

		// opens BAM file
		public void Open(string filename)
		{
			if (!File.Exists(filename))
			{
				throw new ApplicationException(string.Format("ERROR: The supplied BAM filename ({0}) does not exist.",
															 filename));
			}

			BamPath = filename;

			// sanity check: make sure this is a GZIP file
			ushort gzipMagicNumber;
			using (BinaryReader checkReader = new BinaryReader(new FileStream(filename, FileMode.Open, FileAccess.Read, FileShare.Read)))
			{
				gzipMagicNumber = checkReader.ReadUInt16();
			}

			if (gzipMagicNumber != BamConstants.GzipMagicNumber)
			{
				throw new ApplicationException(string.Format(
					"The input file ({0}) does not seem to be a BAM file. A GZIP magic number could not be found.", filename));
			}

			// open our file streams
			try
			{
				BgzfFileStream = new FileStream(filename, FileMode.Open, FileAccess.Read, FileShare.Read);
				_reader = new BinaryReader(BgzfFileStream);
			}
			catch (IOException e)
			{
				throw new ApplicationException(string.Format(
					"ERROR: Unable to open the BAM file ({0}) for reading: {1}", filename, e.Message));
			}

			IsOpen = true;

			// open the index file if it exists
			string indexPath = string.Format("{0}.bai", filename);
			if (_index.ReadIndex(indexPath)) _hasIndex = true;

			LoadHeaderData();
			LoadReferenceData();

			// store file offset of first alignment
			_alignmentsOffset = ((BlockAddress << 16) | ((long)BlockOffset & 0xFFFF));
		}

		// reads data from the BGZF block
		public int Read(ref byte[] data, uint dataLength)
		{
			if (dataLength == 0) return 0;
			if (dataLength > data.Length)
			{
				// 2014.05.14 YH: for very long reads, re-allocate this buffer to twice the data length
				data = new byte[2 * dataLength];
			}
			int outputIndex = 0;
			int numBytesRead = 0;

			while (numBytesRead < dataLength)
			{
				int bytesAvailable = BlockLength - BlockOffset;

				if (bytesAvailable <= 0)
				{
					if (ReadBlock() != 0) return -1;
					bytesAvailable = BlockLength - BlockOffset;
					if (bytesAvailable <= 0) break;
				}

				int copyLength = Math.Min((int)dataLength - numBytesRead, bytesAvailable);
				Buffer.BlockCopy(UncompressedBlock, BlockOffset, data, outputIndex, copyLength);

				BlockOffset += copyLength;
				outputIndex += copyLength;
				numBytesRead += copyLength;
			}

			if (BlockOffset == BlockLength)
			{
				BlockAddress = BgzfFileStream.Position;
				BlockOffset = 0;
				BlockLength = 0;
			}

			return numBytesRead;
		}

		// reads another BGZF block
		public int ReadBlock()
		{
			byte[] header = new byte[BamConstants.BlockHeaderLength];
			long blockAddress = BgzfFileStream.Position;

			int count = _reader.Read(header, 0, BamConstants.BlockHeaderLength);
			if (count == 0)
			{
				BlockLength = 0;
				return 0;
			}

			if (count != BamConstants.BlockHeaderLength)
			{
				throw new ApplicationException(
					string.Format("ERROR: Expected to read {0} bytes from the block header, but only read {1} bytes.",
								  BamConstants.BlockHeaderLength, count));
			}

			int blockLength = BitConverter.ToUInt16(header, 16) + 1;
			int remaining = blockLength - BamConstants.BlockHeaderLength;

			Buffer.BlockCopy(header, 0, CompressedBlock, 0, BamConstants.BlockHeaderLength);
			count = _reader.Read(CompressedBlock, BamConstants.BlockHeaderLength, remaining);

			if (count != remaining)
			{
				throw new ApplicationException(
					string.Format("ERROR: Expected to read {0} bytes from the block header, but only read {1} bytes.",
					remaining, count));
			}

			count = SafeNativeMethods.UncompressBlock(CompressedBlock, (uint)blockLength, UncompressedBlock, MaxBlockSize);
			if (count < 0) return -1;

			if (BlockLength != 0) BlockOffset = 0;

			BlockAddress = blockAddress;
			BlockLength = count;

			return 0;
		}

		// returns file pointer to beginning of alignments
		public void Rewind()
		{
			BlockLength = 0;
			BlockOffset = (int)(_alignmentsOffset & 0xFFFF);
			BlockAddress = (_alignmentsOffset >> 16) & 0xFFFFFFFFFFFF;

			try
			{
				BgzfFileStream.Seek(BlockAddress, SeekOrigin.Begin);
			}
			catch (IOException e)
			{
				throw new ApplicationException(string.Format("ERROR: Unable to seek in the BAM file: {0}", e.Message));
			}
		}

		/// <summary>
		///     sets the file pointer to the specified location (minus some minor bit-shifting)
		/// </summary>
		public void Seek(ulong offset)
		{
			int blockOffset = (int)(offset & 0xFFFF);
			long blockAddress = (long)(offset >> 16);

			BgzfFileStream.Position = blockAddress;

			BlockLength = 0;
			BlockAddress = blockAddress;
			BlockOffset = blockOffset;
		}
	}
}