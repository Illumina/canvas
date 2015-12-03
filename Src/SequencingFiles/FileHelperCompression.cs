using System;
using System.Collections.Generic;
using System.IO;

namespace SequencingFiles
{
    public static class FileHelperCompression
    {

        private static readonly object IndexLock = new object();

        public delegate void ErrorHandler(string message);
        public static event ErrorHandler Error;

        /// <summary>
        /// Function to load a single tile from a lane aggregate, block compressed bcl file 
        /// </summary>
        /// <param name="filePath">Path to bgzf compressed aggregate bcl file</param>
        /// <param name="uncompressedStart">The cluster offset for this tile from bci file</param>
        /// <param name="clusterCount">The number of clusters to read for this tile</param>
        /// <param name="uncompressedBlockSizeByAddress">Optional: If bgzi index table is known, supply it here (set to null otherwise)</param>
        /// <returns></returns>
        public static byte[][] LoadIndexedBgzippedCallsAndScores(string filePath, long uncompressedStart,
                                                                 int clusterCount, Dictionary<long,int> uncompressedBlockSizeByAddress)
        {
            
            // try to read or index the file
            if (uncompressedBlockSizeByAddress == null || uncompressedBlockSizeByAddress.Count == 0)
            {
                uncompressedBlockSizeByAddress = IndexBgzfBlocks(filePath);
            }
            ulong endPosition;
            ulong offset = (ulong)BgzfPositionFromUncompressedPosition(uncompressedStart, uncompressedBlockSizeByAddress);
            return LoadBgzippedCallsAndScores(filePath,offset,clusterCount, out endPosition);

        }

        /// <summary>
        /// Function to index block compressed bgzf file 
        /// Allows indexing into aggregate, block compressed files 
        /// </summary>
        /// <param name="bgzfPath">Path to bgzf file</param>
        /// <returns>Uncompressed block size by seek offset in compressed BCL file</returns>
        public static Dictionary<long, int> IndexBgzfBlocks(string bgzfPath)
        {
            string bgzfIndexPath = bgzfPath.Replace("bgzf", "bgzi");
            Dictionary<long, int> uncompressedBlockSizeByAddress = new Dictionary<long, int>();

            // make sure this hasn't been indexed already
            if (File.Exists(bgzfPath))
            {
                lock (IndexLock) // only one thread at a time should index an aggregate file, and other threads must wait to read until indexing is done
                {
                    if (!File.Exists(bgzfIndexPath))
                    {
                        // note that creating the reader will read the first block in the file
                        BgzipReader reader = new BgzipReader(bgzfPath);

                        while (reader.BlockLength > 0)
                        {
                            uncompressedBlockSizeByAddress.Add(reader.BlockAddress, reader.BlockLength);
                            reader.ReadBlock();
                        }
                        WriteBgzfIndexFile(bgzfIndexPath, uncompressedBlockSizeByAddress);
                    }
                    else
                    {
                        uncompressedBlockSizeByAddress = ReadBgzfIndexFile(bgzfIndexPath);
                    }
                }
            }

            return uncompressedBlockSizeByAddress;
        }

        #region PrivateFunctions

        private static void OnError(string message)
        {
            if (Error != null) Error(message);
        }

        private static byte[][] LoadBgzippedCallsAndScores(string filePath, ulong startPosition, int clusterCount, out ulong endPosition)
        {
            byte[][] callsAndScores;
            using (BgzipReader reader = new BgzipReader(filePath))
            {

                reader.Seek(startPosition+4);
                byte[] buffer = reader.ReadBytes(clusterCount, out endPosition);

                callsAndScores = new byte[2][];
                callsAndScores[0] = new byte[clusterCount];
                callsAndScores[1] = new byte[clusterCount];

                for (int clusterIndex = 0; clusterIndex < clusterCount; clusterIndex++)
                {
                    if (buffer[clusterIndex] != 0) // can't be uncalled, so quality score must be > 0
                        callsAndScores[0][clusterIndex] = (byte) ((buffer[clusterIndex] & 3) + 8);
                    else
                        callsAndScores[0][clusterIndex] = 0; // no-call reserved
                    callsAndScores[1][clusterIndex] = (byte) (buffer[clusterIndex] >> 2);
                }
            }
            return callsAndScores;
        }

        private static Dictionary<long, int> ReadBgzfIndexFile(string bgzfIndexPath)
        {
            Dictionary<long, int> returnDictionary = new Dictionary<long, int>();

            try
            {
                using (FileStream fs = new FileStream(bgzfIndexPath, FileMode.Open, FileAccess.Read,
                                                                FileShare.Read))
                using (BinaryReader br = new BinaryReader(fs))
                {
                    br.ReadByte(); // version
                    int numPairs = br.ReadInt32();

                    for (int pairIndex = 0; pairIndex < numPairs; pairIndex++)
                    {
                        long address = br.ReadInt64();
                        int blockSize = br.ReadInt32();
                        returnDictionary.Add(address,blockSize);
                    }
                }
            }
            catch (Exception exception)
            {
                OnError("FileCompressionHelper Error: " + exception.Message + " StackTrace: " + exception.StackTrace);
            }

            return returnDictionary;
        }

        private static void WriteBgzfIndexFile(string bgzfIndexPath, Dictionary<long, int> uncompressedBlockSizeByAddress)
        {
            try
            {
                using (FileStream fs = new FileStream(bgzfIndexPath, FileMode.Create, FileAccess.Write, FileShare.None))
                using (BinaryWriter bw = new BinaryWriter(fs))
                {
                    // version
                    bw.Write((byte)1);
                    // number of pairs in file
                    bw.Write(uncompressedBlockSizeByAddress.Count);
                    foreach (long address in uncompressedBlockSizeByAddress.Keys)
                    {
                        bw.Write(address);
                        bw.Write(uncompressedBlockSizeByAddress[address]);
                    }
                }
            }
            catch (Exception exception)
            {
                OnError("FileCompressionHelper Error: " + exception.Message + " StackTrace: " + exception.StackTrace);
            }
        }

        private static long BgzfPositionFromUncompressedPosition(long uncompressedPosition, Dictionary<long, int> uncompressedBlockSizeByAddress)
        {
            long returnPosition = 0;

            long uncompressedBytesSoFar = 0;
            foreach (long compressedAddress in uncompressedBlockSizeByAddress.Keys)
            {
                int currentBlockSize = uncompressedBlockSizeByAddress[compressedAddress];
                if (uncompressedBytesSoFar + currentBlockSize > uncompressedPosition)
                {
                    int blockOffset = (int)(uncompressedPosition - uncompressedBytesSoFar);
                    returnPosition = PositionInBgzf(compressedAddress, blockOffset);
                    break;
                }
                uncompressedBytesSoFar += currentBlockSize;
            }

            return returnPosition;

        }

        // use this construction to "Seek" to a position in a bgzf file
        private static long PositionInBgzf(long compressedBlockAddress, int uncompressedBlockOffset)
        {
            return ((compressedBlockAddress << 16) | (long)(uncompressedBlockOffset) & 0xFFFF);
        }

        #endregion
    }
}
