using System;
using System.IO;
using System.Runtime.InteropServices;

namespace Illumina.Zlib
{
    // native calls for a gzfile
    public abstract class ZBase : Stream
    {
        protected static readonly int ZFinish = 4;
        protected static readonly int ZSyncFlush = 2;
        protected static readonly int ZErrno = -1;
        protected static readonly int ZStreamError = -2;
        protected static readonly int ZOk = 0;
        public static readonly string CompressionLevelHuffmanOnly = "wb1h";
        protected IntPtr Gzfile = IntPtr.Zero;
        protected Int64 ProtectedPosition = 0;

        [DllImport("FileCompression.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern int CompressBlock(byte[] uncompressedBlock, ref int uncompressedLen, byte[] compressedBlock,
                                               uint compressedLen, int compressionLevel);

        [DllImport("FileCompression.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern int UncompressBlock(byte[] compressedBlock, uint compressedLen, byte[] uncompressedBlock,
                                                 uint uncompressedLen);

        [DllImport("FileCompression.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern int gzclose(IntPtr fh);

        [DllImport("FileCompression.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern int gzeof(IntPtr fh);

        [DllImport("FileCompression.dll", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi,
            BestFitMapping = false)]
        public static extern IntPtr gzopen(string path, string mode);

        [DllImport("FileCompression.dll", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
        public static extern int gzreadOffset(IntPtr fh, byte[] buffer, int offset, uint len);

        [DllImport("FileCompression.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern int gzrewind(IntPtr fh);

        [DllImport("FileCompression.dll", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
        public static extern int gzwriteOffset(IntPtr fh, byte[] buffer, int offset, uint len);

        [DllImport("FileCompression.dll", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
        public static extern unsafe int gzwriteOffset(IntPtr fh, ushort* buffer, int offset, uint len);

        protected void Error(string reason)
        {
            //Int32 errnum = 0;
            //IntPtr msg = gzerror(gzfile, out errnum);
            //String str = reason;
            //if (errnum != 0)
            //    str = str + " num=" + errnum;

            //if (msg != IntPtr.Zero)
            //    str = str + " desc=" + Marshal.PtrToStringAnsi(msg);

            throw new ApplicationException(reason);
        }
    }

    public class ZInputStream : ZBase
    {
        public ZInputStream(String fname)
        {
            Gzfile = gzopen(fname, "rb");
            if (Gzfile == IntPtr.Zero)
                Error("gzopen fail to read " + fname);
        }

        public override Boolean CanRead
        {
            get { return true; }
        }

        public override Boolean CanSeek
        {
            get { return false; }
        }

        public override Boolean CanWrite
        {
            get { return false; }
        }

        public override Int64 Length
        {
            get { return ProtectedPosition; }
        }

        public override Int64 Position
        {
            get { return ProtectedPosition; }

            set { }
        }

        public override int Read(byte[] bytes, int off, int len)
        {
            int res = gzreadOffset(Gzfile, bytes, off, (uint) len);
            if (res == -1)
                Error("gzread fail");
            ProtectedPosition += res;
            return res;
        }

        public override void Close()
        {
            if (Gzfile == IntPtr.Zero) return;
            int res = gzclose(Gzfile);
            Gzfile = IntPtr.Zero;
            if (res == ZErrno)
                Error("gzclose fail Z_ERRNO");
            if (res == ZStreamError)
                Error("gzclose fail Z_STREAM_ERROR");
            if (res != ZOk)
                Error("gzclose fail unknown");
        }

        public override void Flush()
        {
        }

        public override void Write(Byte[] buffer, Int32 offset, Int32 count)
        {
            throw new NotImplementedException();
        }

        public override void SetLength(Int64 value)
        {
            throw new NotImplementedException();
        }

        public override Int64 Seek(Int64 offset, SeekOrigin origin)
        {
            throw new NotImplementedException();
        }
    }

    public class ZOutputStream : ZBase
    {
        public ZOutputStream(string fname)
        {
            Gzfile = gzopen(fname, "wb"); // huffman only compression, level 1
            if (Gzfile == IntPtr.Zero)
                Error("gzopen write " + fname);
        }

        public ZOutputStream(string fname, string compressionLevel)
        {
            Gzfile = gzopen(fname, compressionLevel); // "wb1h" for huffman only compression, level 1
            if (Gzfile == IntPtr.Zero)
                Error("gzopen write " + fname);
        }

        public override Boolean CanRead
        {
            get { return false; }
        }

        public override Boolean CanSeek
        {
            get { return false; }
        }

        public override Boolean CanWrite
        {
            get { return true; }
        }

        public override Int64 Length
        {
            get { return ProtectedPosition; }
        }

        public override Int64 Position
        {
            get { return ProtectedPosition; }

            set { }
        }

        public override void WriteByte(byte b)
        {
            byte[] b1 = new byte[1];
            b1[0] = b;
            Write(b1, 0, 1);
        }

        public void Write(byte[] b)
        {
            Write(b, 0, b.Length);
        }

        public override void Write(byte[] b, int off, int len)
        {
            int res = gzwriteOffset(Gzfile, b, off, (uint) len);
            if (res == -1)
                Error("gzwrite fail -1");
            if (res != len)
                Error("gzwrite fail partial write");

            ProtectedPosition += res;
        }

        public unsafe void Write(ushort* buffer, int off, int len)
        {
            int res = gzwriteOffset(Gzfile, buffer, off, (uint) len);
            if (res == -1)
                Error("gzwrite fail -1");
            if (res != len)
                Error("gzwrite fail partial write");

            ProtectedPosition += res;
        }

        public override void Close()
        {
            if (Gzfile == IntPtr.Zero) return;
            //gzflush(gzfile, Z_FINISH);
            int res = gzclose(Gzfile);
            Gzfile = IntPtr.Zero;
            if (res == ZErrno)
                Error("gzclose fail Z_ERRNO");
            if (res == ZStreamError)
                Error("gzclose fail Z_STREAM_ERROR");
            if (res != ZOk)
                Error("gzclose fail unknown");
        }

        public override void Flush()
        {
            // No-op, for now:
            //int res = gzflush(gzfile, Z_SYNC_FLUSH);
        }

        public override Int32 Read(Byte[] buffer, Int32 offset, Int32 count)
        {
            throw new NotImplementedException();
        }

        public override void SetLength(Int64 value)
        {
            throw new NotImplementedException();
        }

        public override Int64 Seek(Int64 offset, SeekOrigin origin)
        {
            throw new NotImplementedException();
        }
    }
}