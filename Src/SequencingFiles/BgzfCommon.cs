using System;
using System.IO;
using SequencingFiles.Compression;

namespace SequencingFiles
{
    public abstract class BgzfCommon : ClosableDisposable
    {
        // constants
        protected const int MaxBlockSize = 64 * 1024;

        // variables
        protected internal Stream BgzfFileStream;
        protected internal long BlockAddress;

        protected internal int BlockLength;
        protected internal int BlockOffset;
        public byte[] CompressedBlock;
        protected bool IsOpen = false;

        protected byte[] UncompressedBlock;

        // constructor
        protected BgzfCommon()
        {
            // initialize our buffers
            try
            {
                CompressedBlock = new byte[MaxBlockSize];
                UncompressedBlock = new byte[MaxBlockSize];
            }
            catch (OutOfMemoryException e)
            {
                throw new ApplicationException(
                    string.Format("ERROR: Unable to allocate memory for the compressed and uncompressed blocks: {0}",
                                  e.Message));
            }
        }

        // returns the location of the file pointer
        public ulong Tell()
        {
            return (((ulong)BlockAddress << 16) | ((ulong)BlockOffset & 0xFFFF));
        }
    }

    /// <summary>
    ///     Common base for bgzf file writers (bam and bgzip)
    /// </summary>
    public class BgzfWriterCommon : BgzfCommon
    {
        #region member variables

        private readonly int _compressionLevel;
        private BinaryWriter _writer;

        #endregion

        public BgzfWriterCommon(int compressionLevel)
        {
            _compressionLevel = compressionLevel;
        }

        protected void Open(string filename)
        {
            if (IsOpen) Close();

            BgzfFileStream = new FileStream(filename, FileMode.Create);
            _writer = new BinaryWriter(BgzfFileStream);
            IsOpen = true;
        }

        protected void Open(Stream outStream)
        {
            if (IsOpen) Close();

            BgzfFileStream = outStream;
            _writer = new BinaryWriter(BgzfFileStream);
            IsOpen = true;
        }

        // close Bgzipped file
        public override void Close()
        {
            if (!IsOpen) return;
            IsOpen = false;

            // flush the current BGZF block
            FlushBlock();

            // write an empty block (as EOF marker)
            FlushSingleBlock();

            _writer.Close();
            BgzfFileStream.Close();
        }

        private int FlushSingleBlock()
        {
            int blockLength = SafeNativeMethods.CompressBlock(
                UncompressedBlock,
                ref BlockOffset,
                CompressedBlock,
                MaxBlockSize,
                _compressionLevel);

            try
            {
                _writer.Write(CompressedBlock, 0, blockLength);
            }
            catch (IOException e)
            {
                throw new ApplicationException(
                    string.Format("ERROR: An IO exception occurred when flushing the BGZF block: {0}", e.Message));
            }

            return blockLength;
        }

        // flushes the data in the BGZF block
        protected void FlushBlock()
        {
            // flush all of the remaining blocks
            while (BlockOffset > 0)
            {
                BlockAddress += FlushSingleBlock();
            }
        }

        // writes data to the BGZF block
        protected int Write(byte[] data, uint dataLength)
        {
            // initialize
            const int blockLength = MaxBlockSize;
            int numBytesWritten = 0;
            int inputIndex = 0;

            // copy the data to the buffer
            while (numBytesWritten < dataLength)
            {
                int copyLength = Math.Min(blockLength - BlockOffset, (int)dataLength - numBytesWritten);

                Buffer.BlockCopy(data, inputIndex, UncompressedBlock, BlockOffset, copyLength);

                BlockOffset += copyLength;
                inputIndex += copyLength;
                numBytesWritten += copyLength;

                if (BlockOffset == blockLength) FlushBlock();
            }

            return numBytesWritten;
        }
    }
}