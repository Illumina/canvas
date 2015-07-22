using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using SequencingFiles.Compression;


namespace SequencingFiles
{
    /// <summary>
    ///     Reads from bgzip format
    ///     Any valid bgzip file can be read with a gzip reader but this class brings vcf indexing into C#
    ///     through the Position() and Seek() methods.
    /// </summary>
    public class BgzipReader : BgzfCommon
    {
        #region member variables 
        protected const byte GCLineFeed = 10;
        internal protected byte[] LineBuffer = new byte[MaxBlockSize * 2];
        private byte[] Line = new byte[MaxBlockSize * 2];
        private BinaryReader _reader;
        #endregion

        public BgzipReader(string filename)
            : base()
        {
            //Console.WriteLine(filename);
            Open(filename);
        }

        private void Open(string filename)
        {
            if (IsOpen) Close();

            BgzfFileStream = new FileStream(filename, FileMode.Open, FileAccess.Read, FileShare.Read);
            _reader = new BinaryReader(BgzfFileStream);
            IsOpen = true;
            FillBuffer();
        }

        // close Bgzipped file
        public override void Close()
        {
            if (!IsOpen) return;
            IsOpen = false;

            _reader.Close();
            BgzfFileStream.Close();
        }

        /// <summary>
        ///     Gets the next string in the file.
        /// </summary>
        /// <returns>Returns the next string or null.</returns>
        public string ReadLine()
        {
            string s = "";

            // skip if the file is not currently open or if we don't have any data in the buffer
            if (!IsOpen || (BlockLength <= 0)) return null;

            int crOffset = -1;
            while (true)
            {
                crOffset = Array.IndexOf(LineBuffer, GCLineFeed, BlockOffset, BlockLength - BlockOffset);
                if (crOffset != -1)
                {
                    s = s + Encoding.ASCII.GetString(LineBuffer, BlockOffset, crOffset - BlockOffset);
                    BlockOffset = crOffset + 1;
                    break;
                }
                else
                {
                    int remainingLen = BlockLength - BlockOffset;

                    s = s + Encoding.ASCII.GetString(LineBuffer, BlockOffset, remainingLen);

                    FillBuffer();

                    if (BlockLength <= 0)
                    {
                        //end of file?
                        return string.IsNullOrEmpty(s) ? null : s;
                    }
                }
            }
            return s;
        }

        /// <summary>
        ///     Fills the buffer
        /// </summary>
        /// <param name="bufferOffset"></param>
        public void FillBuffer()
        {
            if (ReadBlock() != 0) return;
            Buffer.BlockCopy(UncompressedBlock, 0, LineBuffer, 0, BlockLength);
            BlockOffset = 0;
            return;
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

            // total block size - 1 is located in header starting at byte 16 (zero indexed)
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

            // if BlockLength is zero we must have done a seek so don't update BlockOffset
            if (BlockLength != 0) BlockOffset = 0;

            BlockAddress = blockAddress;
            BlockLength = count;
            return 0;
        }

        /// <summary>
        /// Read and return byteCount bytes, starting from the current position, and returning the new position.
        /// </summary>
        public byte[] ReadBytes(int byteCount, out ulong nextPosition)
        {
            byte[] outputBuffer = new byte[byteCount];
            int outputBufferPos = 0;
            while (outputBufferPos < outputBuffer.Length)
            {
                // Add bytes from the current block, if available:
                if (this.BlockOffset < this.BlockLength)
                {
                    int copyBytes = Math.Min(outputBuffer.Length - outputBufferPos, this.BlockLength - this.BlockOffset);
                    Array.Copy(this.UncompressedBlock, this.BlockOffset, outputBuffer, outputBufferPos, copyBytes);
                    this.BlockOffset += copyBytes;
                    outputBufferPos += copyBytes;
                }
                else
                {
                    int readBytes = this.ReadBlock();
                    if (BlockLength == 0)
                    {
                        throw new Exception(string.Format("Error: Attempted to read beyond the end of a bgzipped file; file may be truncated",
                            this.Position()));
                    }
                }
            }
            nextPosition = (ulong)this.Position();
            return outputBuffer;
        }

        public long Position()
        {
            return (BlockAddress) << 16 | ((long)(BlockOffset) & 0xFFFF);
        }

        /// <summary>
        ///     sets the file pointer to the specified location (minus some minor bit-shifting)
        /// </summary>
        public void Seek(ulong offset)
        {
            BlockAddress = (long)(offset >> 16);
            BgzfFileStream.Position = BlockAddress;
            FillBuffer();
            BlockOffset = (int)(offset & 0xFFFF);
        }
    }
}