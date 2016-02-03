using System;
using System.IO;
using SequencingFiles.Compression;
using System.Collections.Concurrent;
using System.Threading;
using System.Collections.Generic;

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

        // Stores the compressed data. CompressedData objects are queued up to be
        // written to file on a single thread.
        private class CompressedData
        {
            public CompressedData()
            {
                _buffers = new List<byte[]>();
                _bufferSizes = new List<int>();
                _autoEvent = new AutoResetEvent(false);
            }
            public List<byte[]> _buffers;
            public List<int> _bufferSizes;

            // Set this flag when the compressed data is ready to be written.
            // This synchronizes the writing thread and processing threads.
            public AutoResetEvent _autoEvent;
        }

        // Stores the uncompressed data. Data objects are queued up to be processed
        // on separate threads.
        private class Data
        {
            public Data(byte[] buffer, int bufferSize, CompressedData compressedData)
            {
                _buffer = buffer;
                _bufferSize = bufferSize;
                _compressedData = compressedData;
            }
            public byte[] _buffer;
            public int _bufferSize;

            public CompressedData _compressedData;
        }

        protected readonly int CompressionLevel;
        protected readonly int CompressionStrategy;

        // The compression threads pull the data off of this queue.
        private BlockingCollection<Data> _compressQueue;

        // The writing thread pull compressed blocks off of this queue.
        // Note that the order blocks are compressed in is nondeterministic,
        // but the writing is always done in fixed order by using the _autoEvent signal.
        private BlockingCollection<CompressedData> _writeQueue;

        private Thread[] _compressionThreads;
        private Thread _writingThread;
        private BinaryWriter _writer;
        protected int NumThreads;

        #endregion

        public BgzfWriterCommon(int compressionLevel, int compressionStrategy, int numThreads = 1)
        {
            CompressionLevel = compressionLevel;
            CompressionStrategy = compressionStrategy;
            NumThreads = numThreads;

            _compressQueue = null;
            _writeQueue = null;
            _compressionThreads = null;
            _writingThread = null;
        }

        private void InitThreads()
        {
            if (NumThreads > 1)
            {
                _compressQueue = new BlockingCollection<Data>(2 * NumThreads);
                _writeQueue = new BlockingCollection<CompressedData>(2 * NumThreads);
                _compressionThreads = new Thread[NumThreads];

                for (int i = 0; i < NumThreads; ++i)
                {
                    _compressionThreads[i] = new Thread(() => CompressionThread());
                    _compressionThreads[i].Name = string.Format("Compression thread {0}", i);
                    _compressionThreads[i].Start();
                }

                _writingThread = new Thread(() => WritingThread());
                _writingThread.Name = string.Format("Writing thread");
                _writingThread.Start();
            }
        }

        private void CompressionThread()
        {
            try
            {
                while (true)
                {
                    Data uncompressedBlock = _compressQueue.Take();
                    int blockOffset = uncompressedBlock._bufferSize;
                    while (blockOffset > 0)
                    {
                        byte[] compressedBlock = new byte[MaxBlockSize];
                        int blockLength = SafeNativeMethods.CompressBlock(
                            uncompressedBlock._buffer,
                            ref blockOffset,
                            compressedBlock,
                            MaxBlockSize,
                            CompressionLevel,
                            CompressionStrategy);

                        uncompressedBlock._compressedData._buffers.Add(compressedBlock);
                        uncompressedBlock._compressedData._bufferSizes.Add(blockLength);
                    }

                    // Finished compressing the block. Signal that it's ok to write it to disk.
                    uncompressedBlock._compressedData._autoEvent.Set();
                }
            }
            catch (InvalidOperationException)
            {
                // We're done. There are no more blocks to process.
            }
        }

        private void WritingThread()
        {
            try
            {
                while (true)
                {
                    CompressedData compressedBlock = _writeQueue.Take();

                    // Wait for the block to be compressed.
                    compressedBlock._autoEvent.WaitOne();
                    for (int i = 0; i < compressedBlock._buffers.Count; ++i)
                    {
                        _writer.Write(compressedBlock._buffers[i], 0, compressedBlock._bufferSizes[i]);
                    }
                }
            }
            catch (IOException e)
            {
                throw new ApplicationException(
                    string.Format("ERROR: An IO exception occurred when flushing the BGZF block: {0}", e.Message));
            }
            catch (InvalidOperationException)
            {
                // We're done. There are no more blocks to process.
            }
        }

        protected void Open(string filename)
        {
            if (IsOpen) Close();

            BgzfFileStream = new FileStream(filename, FileMode.Create);
            _writer = new BinaryWriter(BgzfFileStream);
            IsOpen = true;
            InitThreads();
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

            if (_compressQueue != null)
            {
                // Signal that we're not going to add any more data to the queue.
                _compressQueue.CompleteAdding();

                // Wait for our compression threads to finish.
                int numCompressionThreads = _compressionThreads.Length;
                for (int i = 0; i < numCompressionThreads; ++i)
                {
                    _compressionThreads[i].Join();
                }

                // Signal to the writing thread that compression is complete.
                _writeQueue.CompleteAdding();

                // Wait for the writing thread to finish.
                _writingThread.Join();
            }

            // write an empty block (as EOF marker)
            FlushSingleBlock();

            _writer.Close();
            BgzfFileStream.Close();

            _compressQueue = null;
            _compressionThreads = null;
            _writeQueue = null;
            _writingThread = null;
        }

        private int FlushSingleBlock()
        {
            int blockLength = SafeNativeMethods.CompressBlock(UncompressedBlock,
                ref BlockOffset, CompressedBlock, MaxBlockSize,
                CompressionLevel, CompressionStrategy);

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
            if (_compressQueue != null)
            {
                // Save the block to the queue to be processed on multiple threads.
                CompressedData compressedData = new CompressedData();
                _writeQueue.Add(compressedData);
                _compressQueue.Add(new Data((byte[])UncompressedBlock.Clone(), BlockOffset, compressedData));
                BlockOffset = 0;
            }
            else
            {
                // flush all of the remaining blocks
                while (BlockOffset > 0)
                {
                    BlockAddress += FlushSingleBlock();
                }
            }
        }

        // writes data to the BGZF block
        public int Write(byte[] data, uint dataLength)
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

                if (BlockOffset == blockLength)
                {
                    FlushBlock();
                }
            }

            return numBytesWritten;
        }
    }
}