using System;
using System.Collections.Generic;
using System.Text;
using System.IO;

namespace SequencingFiles
{
    public class BamWriter : BgzfWriterCommon
    {
        private static class ConversionHelper
        {
            static public uint[] BaseToNumber;
            static public uint[] CigarOpToNumber;

            static ConversionHelper()
            {
                // create our base LUT
                BaseToNumber = new uint[256];
                for (int i = 0; i < 256; i++) BaseToNumber[i] = BamConstants.LutError;

                BaseToNumber['='] = 0;
                BaseToNumber['a'] = 1;
                BaseToNumber['A'] = 1;
                BaseToNumber['c'] = 2;
                BaseToNumber['C'] = 2;
                BaseToNumber['g'] = 4;
                BaseToNumber['G'] = 4;
                BaseToNumber['t'] = 8;
                BaseToNumber['T'] = 8;
                BaseToNumber['n'] = 15;
                BaseToNumber['N'] = 15;

                // create our CIGAR LUT
                CigarOpToNumber = new uint[256];
                for (int i = 0; i < 256; i++) CigarOpToNumber[i] = BamConstants.LutError;
                CigarOpToNumber['m'] = 0;
                CigarOpToNumber['M'] = 0;
                CigarOpToNumber['i'] = 1;
                CigarOpToNumber['I'] = 1;
                CigarOpToNumber['d'] = 2;
                CigarOpToNumber['D'] = 2;
                CigarOpToNumber['n'] = 3;
                CigarOpToNumber['N'] = 3;
                CigarOpToNumber['s'] = 4;
                CigarOpToNumber['S'] = 4;
                CigarOpToNumber['h'] = 5;
                CigarOpToNumber['H'] = 5;
                CigarOpToNumber['p'] = 6;
                CigarOpToNumber['P'] = 6;
                CigarOpToNumber['='] = 7;
                CigarOpToNumber['x'] = 8;
                CigarOpToNumber['X'] = 8;
            }
        }
        #region member variables

        private byte[] _outputBuffer;
        protected string SamHeader;
        protected List<GenomeMetadata.SequenceMetadata> References;
        protected string Filename;

        #endregion

        // constructor
        public BamWriter(int compressionLevel = BamConstants.DefaultCompression, int numThreads = 1)
            : this(null, null, null, compressionLevel, numThreads)
        {
            _outputBuffer = null;
        }

        // constructor
        public BamWriter(string filename, string samHeader, List<GenomeMetadata.SequenceMetadata> references, int compressionLevel = BamConstants.DefaultCompression, int numThreads = 1)
            : base(compressionLevel, numThreads)
        {
            Filename = filename;
            SamHeader = samHeader;
            References = references;
            if (filename != null) Open(filename, samHeader, references);
        }

        private void Initialize()
        {
            if (_outputBuffer == null)
            {
                _outputBuffer = new byte[4096];
            }
        }

        // opens BAM file
        public void Open(string filename, string samHeader, List<GenomeMetadata.SequenceMetadata> references)
        {
            Initialize();
            base.Open(filename);

            // ================
            // write the header
            // ================

            // write the BAM signature
            byte[] buffer = Encoding.ASCII.GetBytes(BamConstants.MagicNumber);
            Write(buffer, (uint)BamConstants.MagicNumber.Length);

            // write the SAM header
            uint samHeaderLen = (uint)samHeader.Length;
            buffer = BitConverter.GetBytes(samHeaderLen);
            Write(buffer, 4);

            if (!string.IsNullOrEmpty(samHeader))
            {
                buffer = Encoding.ASCII.GetBytes(samHeader);
                Write(buffer, samHeaderLen);
            }

            // write the number of reference sequences
            uint numReferenceSequences = (uint)references.Count;
            buffer = BitConverter.GetBytes(numReferenceSequences);
            Write(buffer, 4);

            // =============================
            // write the sequence dictionary
            // =============================

            foreach (GenomeMetadata.SequenceMetadata refSeq in references)
            {
                // write the reference sequence name
                string referenceName = refSeq.Name + '\0';
                uint referenceSequenceNameLen = (uint)referenceName.Length;
                buffer = BitConverter.GetBytes(referenceSequenceNameLen);
                Write(buffer, 4);

                buffer = Encoding.ASCII.GetBytes(referenceName);
                Write(buffer, referenceSequenceNameLen);

                // write the reference sequence length
                buffer = BitConverter.GetBytes(refSeq.Length);
                Write(buffer, 4);
            }

            // flush the buffer after the header has been written
            FlushBlock();
        }

        // encodes the supplied query sequence into 4-bit notation
        private void PackBases(ref int offset, uint numEncodedBases, string s)
        {
            for (uint baseIndex = 0; baseIndex < numEncodedBases; baseIndex++) _outputBuffer[offset + baseIndex] = 0;
            byte shift = 4;

            int endOffset = offset + (int)numEncodedBases;

            foreach (char c in s)
            {
                uint baseCode = ConversionHelper.BaseToNumber[c];
                if (baseCode == BamConstants.LutError)
                {
                    throw new ApplicationException(
                        string.Format("ERROR: Encountered an unexpected base ({0}) when packing the read.", c));
                }

                _outputBuffer[offset] |= (byte)(baseCode << shift);
                if (shift == 0) offset++;
                shift ^= 4;
            }

            offset = endOffset;
        }

        // encodes the supplied query sequence into 4-bit notation
        static private void PackBases(ref int offset, ref byte[] buffer, uint numEncodedBases, string s)
        {
            //for (uint baseIndex = 0; baseIndex < numEncodedBases; baseIndex++) buffer[offset + baseIndex] = 0;
            byte shift = 4;

            int endOffset = offset + (int)numEncodedBases;

            foreach (char c in s)
            {
                uint baseCode = ConversionHelper.BaseToNumber[c];
                if (baseCode == BamConstants.LutError)
                {
                    throw new ApplicationException(
                        string.Format("ERROR: Encountered an unexpected base ({0}) when packing the read.", c));
                }

                buffer[offset] |= (byte)(baseCode << shift);
                if (shift == 0) offset++;
                shift ^= 4;
            }

            offset = endOffset;
        }

        static private void PackCigar(ref int offset, ref byte[] buffer, CigarAlignment cigarOps)
        {
            // pack the cigar data into the string

            foreach (CigarOp op in cigarOps)
            {
                uint cigarOp = ConversionHelper.CigarOpToNumber[op.Type];
                if (cigarOp == BamConstants.LutError)
                {
                    throw new ApplicationException(
                        string.Format("ERROR: Encountered an unexpected CIGAR operation ({0}).", op.Type));
                }

                BinaryIO.AddUIntBytes(ref buffer, ref offset, op.Length << BamConstants.CigarShift | cigarOp);
            }
        }


        /// <summary>
        /// Serialize alignment to a byte array, for later flushing to output file.
        /// </summary>
        static public byte[] SerializeAlignment(ref BamAlignment al)
        {
            // initialize
            uint nameLen = (uint)al.Name.Length + 1;
            uint numBases = (uint)al.Bases.Length;
            uint numCigarOperations = (uint)al.CigarData.Count;
            uint packedCigarLen = numCigarOperations * 4;
            uint numEncodedBases = (uint)((numBases / 2.0) + 0.5);
            uint tagDataLen = (uint)al.TagData.Length;
            uint dataBlockSize = nameLen + packedCigarLen + numEncodedBases + numBases + tagDataLen;
            uint alignBlockSize = BamConstants.CoreAlignmentDataLen + dataBlockSize;
            uint blockSize = alignBlockSize + 4;
            byte[] buffer = new byte[blockSize];
            int offset = 0;

            // store the block size
            BinaryIO.AddUIntBytes(ref buffer, ref offset, alignBlockSize);

            // store the BAM core data
            BinaryIO.AddIntBytes(ref buffer, ref offset, al.RefID);
            BinaryIO.AddIntBytes(ref buffer, ref offset, al.Position);
            BinaryIO.AddUIntBytes(ref buffer, ref offset, (al.Bin << 16) | (al.MapQuality << 8) | nameLen);
            BinaryIO.AddUIntBytes(ref buffer, ref offset, (al.AlignmentFlag << 16) | numCigarOperations);
            BinaryIO.AddUIntBytes(ref buffer, ref offset, numBases);
            BinaryIO.AddIntBytes(ref buffer, ref offset, al.MateRefID);
            BinaryIO.AddIntBytes(ref buffer, ref offset, al.MatePosition);
            BinaryIO.AddIntBytes(ref buffer, ref offset, al.FragmentLength);

            // store the alignment name
            BinaryIO.AddNullTerminatedString(ref buffer, ref offset, al.Name);

            // store the packed CIGAR string and packed bases
            PackCigar(ref offset, ref buffer, al.CigarData);
            PackBases(ref offset, ref buffer, numEncodedBases, al.Bases);

            // store the base qualities
            Buffer.BlockCopy(al.Qualities, 0, buffer, offset, al.Qualities.Length);
            offset += al.Qualities.Length;

            // store the tag data
            Buffer.BlockCopy(al.TagData, 0, buffer, offset, al.TagData.Length);
            offset += al.TagData.Length;

            return buffer;
        }

        public void WriteAlignment(byte[] buffer)
        {
            if (!IsOpen)
            {
                throw new ApplicationException(string.Format("ERROR: Tried to write an alignment but the file has not been opened yet."));
            }
            // test if we should flush the block
            if ((BlockOffset + buffer.Length) > MaxBlockSize) FlushBlock();
            Write(buffer, (uint)buffer.Length);
        }

        // writes an alignment
        public void WriteAlignment(BamAlignment al)
        {
            if (!IsOpen)
            {
                throw new ApplicationException(string.Format("ERROR: Tried to write an alignment but the file has not been opened yet."));
            }

            // initialize
            uint nameLen = (uint)al.Name.Length + 1;
            uint numBases = (uint)al.Bases.Length;
            uint numCigarOperations = (uint)al.CigarData.Count;
            uint packedCigarLen = numCigarOperations * 4;
            uint numEncodedBases = (uint)((numBases / 2.0) + 0.5);
            uint tagDataLen = (uint)al.TagData.Length;
            uint dataBlockSize = nameLen + packedCigarLen + numEncodedBases + numBases + tagDataLen;
            uint alignBlockSize = BamConstants.CoreAlignmentDataLen + dataBlockSize;
            uint blockSize = alignBlockSize + 4;
            int offset = 0;

            // test if we should flush the block
            if ((BlockOffset + blockSize) > MaxBlockSize) FlushBlock();

            // redimension the buffer if needed
            if (blockSize > _outputBuffer.Length) _outputBuffer = new byte[blockSize + 1024];

            // store the block size
            BinaryIO.AddUIntBytes(ref _outputBuffer, ref offset, alignBlockSize);

            // store the BAM core data
            BinaryIO.AddIntBytes(ref _outputBuffer, ref offset, al.RefID);
            BinaryIO.AddIntBytes(ref _outputBuffer, ref offset, al.Position);
            BinaryIO.AddUIntBytes(ref _outputBuffer, ref offset, (al.Bin << 16) | (al.MapQuality << 8) | nameLen);
            BinaryIO.AddUIntBytes(ref _outputBuffer, ref offset, (al.AlignmentFlag << 16) | numCigarOperations);
            BinaryIO.AddUIntBytes(ref _outputBuffer, ref offset, numBases);
            BinaryIO.AddIntBytes(ref _outputBuffer, ref offset, al.MateRefID);
            BinaryIO.AddIntBytes(ref _outputBuffer, ref offset, al.MatePosition);
            BinaryIO.AddIntBytes(ref _outputBuffer, ref offset, al.FragmentLength);

            // store the alignment name
            BinaryIO.AddNullTerminatedString(ref _outputBuffer, ref offset, al.Name);

            // store the packed CIGAR string and packed bases
            PackCigar(ref offset, ref _outputBuffer, al.CigarData);
            PackBases(ref offset, numEncodedBases, al.Bases);

            // store the base qualities
            Buffer.BlockCopy(al.Qualities, 0, _outputBuffer, offset, al.Qualities.Length);
            offset += al.Qualities.Length;

            // store the tag data
            Buffer.BlockCopy(al.TagData, 0, _outputBuffer, offset, al.TagData.Length);
            offset += al.TagData.Length;

            // write the alignment
            Write(_outputBuffer, blockSize);
        }
    }
}