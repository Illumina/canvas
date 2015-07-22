using System;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using System.Runtime.InteropServices;

namespace Illumina.Common
{
    /// <summary>
    ///     A point in 2-dimensional space, using floats to store the X and Y coordinates
    /// </summary>
    [StructLayout(LayoutKind.Sequential, Pack = 1)]
    public struct FloatPoint
    {
        public float X;
        public float Y;

        public FloatPoint(float x, float y)
        {
            X = x;
            Y = y;
        }

        /// <summary>
        ///     Load an array of locations (FloatPoint instances) from a binary file.
        ///     The file format includes 8 bytes (unused), then a point count (int32), then each location as a pair
        ///     of floats.
        /// </summary>
        public static FloatPoint[] LoadLocations(string filePath)
        {
            if (!File.Exists(filePath))
            {
                string cFilePath = filePath.Replace("locs", "clocs");

                if (File.Exists(cFilePath))
                    return LoadAndConvertCompressedLocs(cFilePath);

                return new FloatPoint[0];
            }
            FloatPoint[] Locations = new FloatPoint[0];

            using (Stream inputStream = new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                using (BinaryReader reader = new BinaryReader(inputStream))
                {
                    reader.ReadInt32();
                    reader.ReadInt32();
                    int pointCount = reader.ReadInt32();
                    Locations = new FloatPoint[pointCount];
                    byte[] buffer = reader.ReadBytes(pointCount*sizeof (float)*2);
                    float[] tempFloatBuffer = new float[pointCount*2];
                    System.Buffer.BlockCopy(buffer, 0, tempFloatBuffer, 0, buffer.Length);
                    for (int pointIndex = 0; pointIndex < pointCount; pointIndex++)
                    {
                        Locations[pointIndex] = new FloatPoint
                        {
                            X = tempFloatBuffer[pointIndex*2],
                            Y = tempFloatBuffer[pointIndex*2 + 1]
                        };
                    }
                }
            }
            return Locations;
        }

        /// <summary>
        ///     Write out locations to a .locs file.
        /// </summary>
        public static void SaveLocations(FloatPoint[] locs, string outputPath)
        {
            using (FileStream outputStream = new FileStream(outputPath, FileMode.Create))
            {
                using (BinaryWriter locsWriter = new BinaryWriter(outputStream))
                {
                    locsWriter.Write(1);
                    locsWriter.Write(1.0f);
                    if (locs == null)
                    {
                        locsWriter.Write(0);
                    }
                    else
                    {
                        locsWriter.Write(locs.Length);
                        for (int pointIndex = 0; pointIndex < locs.Length; pointIndex++)
                        {
                            locsWriter.Write(locs[pointIndex].X);
                            locsWriter.Write(locs[pointIndex].Y);
                        }
                    }
                }
            }
        }

        private static readonly int blockSize = 25;
        private static readonly int ImageWidth = 2048;
        private static readonly int ImageHeight = 20000;

        /// <summary>
        ///     Load an array of locations (x, y byte pairs) and an array of bin counts from a binary file.
        ///     The file format includes 1 byte (unused), then a number of bins count (uint32), then each bin count followed by
        ///     the locations in the bin as an x, y pair of bytes.
        /// </summary>
        public static byte[] LoadCompressedLocations(string filePath, out byte[] bins)
        {
            if (!File.Exists(filePath))
            {
                bins = new byte[0];
                return new byte[0];
            }

            byte[] content;

            using (FileStream fs = new FileStream(filePath, FileMode.Open, FileAccess.Read))
            {
                content = new byte[fs.Length];
                fs.Read(content, 0, content.Length);
            }

            //byte version = content[0];
            uint numBins = BitConverter.ToUInt32(content, 1);

            int contentOffset = sizeof (byte) + sizeof (uint);

            if (numBins < 0)
                numBins = 0; // just in case data is corrupt

            // get number of locs
            bins = new byte[numBins];
            int numLocs = (content.Length - contentOffset - (int) numBins)/2;

            byte[] locs = new byte[numLocs*2];

            int locsOffset = 0;

            // iterate through and load data
            for (int binIndex = 0; binIndex < numBins; binIndex++)
            {
                bins[binIndex] = content[contentOffset];
                contentOffset += 1;

                for (int index = 0; index < bins[binIndex]; index++)
                {
                    if (contentOffset + (index*2) + 1 < content.Length)
                    {
                        locs[locsOffset] = content[contentOffset + (index*2)];
                        locs[locsOffset + 1] = content[contentOffset + (index*2) + 1];

                        locsOffset += 2;
                    }
                }
                contentOffset += (bins[binIndex]*2);
            }

            return locs;
        }

        /// <summary>
        ///     Write out compressed locations to a .clocs file.
        /// </summary>
        public static void SaveCompressedLocations(byte[] bins, byte[] locs, string outputPath)
        {
            using (FileStream fs = new FileStream(outputPath, FileMode.Create))
            {
                fs.WriteByte(1); // version

                if (locs == null)
                {
                    fs.Write(BitConverter.GetBytes(0), 0, sizeof (ushort));
                }
                else
                {
                    int locsoffset = 0;

                    fs.Write(BitConverter.GetBytes((uint) bins.Length), 0, sizeof (uint));

                    foreach (byte numberPoints in bins)
                    {
                        fs.WriteByte(numberPoints);

                        for (int locsIndex = 0; locsIndex < numberPoints; locsIndex++)
                        {
                            if ((locsIndex*2) + locsoffset + 1 < locs.Length)
                            {
                                fs.WriteByte(locs[(locsIndex*2) + locsoffset]);
                                fs.WriteByte(locs[(locsIndex*2) + locsoffset + 1]);
                            }
                        }

                        locsoffset += (numberPoints*2);
                    }
                }
            }
        }

        public override String ToString()
        {
            return String.Format("({0}, {1}) ", X, Y);
        }

        public static FloatPoint[] DecompressLocs(byte[] bins, byte[] clocs)
        {
            int maxXbins = (int) Math.Ceiling(ImageWidth/(double) blockSize);
            int numLocs = 0;

            // get number of locs
            for (int binIndex = 0; binIndex < bins.Length; binIndex++)
                numLocs += bins[binIndex];

            if (numLocs <= 0)
                return new FloatPoint[0];

            FloatPoint[] locations = new FloatPoint[numLocs];

            float xoffset = 0; // x offset for bin
            float yoffset = 0; // y offset for bin
            int locsoffset = 0;
            int locsIndex = 0;

            for (int binIndex = 0; binIndex < bins.Length; binIndex++)
            {
                ushort numberLocs = Convert.ToUInt16(bins[binIndex]);

                for (int clocsIndex = 0; clocsIndex < numberLocs; clocsIndex++)
                {
                    FloatPoint point = new FloatPoint
                    {
                        X = Convert.ToSingle(clocs[(clocsIndex*2) + locsoffset])/10f + xoffset,
                        Y = Convert.ToSingle(clocs[(clocsIndex*2) + locsoffset + 1])/ 10f + yoffset
                    };

                    locations[locsIndex] = point;

                    locsIndex++;
                }

                locsoffset += (numberLocs*2);

                // update x, y bin offsets for bin
                if (binIndex > 0 && (binIndex + 1)%maxXbins == 0)
                {
                    xoffset = 0;
                    yoffset += blockSize;
                }
                else
                    xoffset += blockSize;
            }

            return locations;
        }

        public static byte[] CompressLocs(FloatPoint[] locs, out byte[] bins, out uint[] origIndexMap)
        {
            int maxXbins = (int) Math.Ceiling(ImageWidth/(double) blockSize);
            int maxYbins = (int) Math.Ceiling(ImageHeight/(double) blockSize);

            bins = new byte[maxXbins*maxYbins];
            List<byte>[] binLocs = new List<byte>[bins.Length];
            List<uint>[] origLocsIndex = new List<uint>[bins.Length];
            int totalLocs = 0;

            // bin locs
            for (int clusterIndex = 0; clusterIndex < locs.Length; clusterIndex++)
            {
                if (locs[clusterIndex].X < 0 || locs[clusterIndex].Y < 0)
                    continue; // out of defined range, can throw away these outliers

                // determine bin index
                int xbin = (int) Math.Floor(locs[clusterIndex].X/blockSize);
                int ybin = (int) Math.Floor(locs[clusterIndex].Y/blockSize);

                int binIndex = ybin*maxXbins + xbin;

                if (xbin >= maxXbins)
                    continue;
                if (ybin >= maxYbins)
                    continue;

                // add locs to bin
                bins[binIndex]++;
                totalLocs++;

                if (binLocs[binIndex] == null)
                    binLocs[binIndex] = new List<byte>(blockSize*blockSize);
                if (origLocsIndex[binIndex] == null)
                    origLocsIndex[binIndex] = new List<uint>(blockSize*blockSize);

                binLocs[binIndex].Add(Convert.ToByte((locs[clusterIndex].X%blockSize)*10f));
                binLocs[binIndex].Add(Convert.ToByte((locs[clusterIndex].Y%blockSize)*10f));
                origLocsIndex[binIndex].Add((uint) clusterIndex);
            }

            // combine all binned locs into one big array
            byte[] clocs = new byte[totalLocs*2];
            origIndexMap = new uint[totalLocs];

            int offset = 0;
            int origOffset = 0;

            for (int binIndex = 0; binIndex < binLocs.Length; binIndex++)
            {
                if (binLocs[binIndex] == null)
                    continue; // no locs in bin, move on

                foreach (byte locByte in binLocs[binIndex])
                {
                    clocs[offset] = locByte;
                    offset++;
                }

                foreach (uint origIndexVal in origLocsIndex[binIndex])
                {
                    origIndexMap[origOffset] = origIndexVal;
                    origOffset++;
                }
            }

            return clocs;
        }

        public static FloatPoint[] LoadAndConvertCompressedLocs(string filepath)
        {
            byte[] bins;
            byte[] clocs = LoadCompressedLocations(filepath, out bins);
            return DecompressLocs(bins, clocs);
        }
    }

    /// <summary>
    ///     A point in 3-dimensional space, using floats to store the X and Y coordinates
    /// </summary>
    [StructLayout(LayoutKind.Sequential, Pack = 1)]
    public struct XYZFloatPoint
    {
        public float X;
        public float Y;
        public float Z;

        public XYZFloatPoint(float x, float y, float z)
        {
            X = x;
            Y = y;
            Z = z;
        }

        public XYZFloatPoint(PointF orig, float z)
        {
            X = orig.X;
            Y = orig.Y;
            Z = z;
        }

        public override int GetHashCode()
        {
            return (X.GetHashCode() ^ Y.GetHashCode() ^ Z.GetHashCode());
        }

        public override bool Equals(object obj)
        {
            XYZFloatPoint other = (XYZFloatPoint) obj;

            return (
                       (X == other.X) &&
                       (Y == other.Y) &&
                       (Z == other.Z)
                   );
        }

        public static bool operator ==(XYZFloatPoint left, XYZFloatPoint right)
        {
            return left.Equals(right);
        }

        public static bool operator !=(XYZFloatPoint left, XYZFloatPoint right)
        {
            return !(left == right);
        }
    }
}