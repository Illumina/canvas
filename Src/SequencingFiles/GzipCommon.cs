using System;
using SequencingFiles.Compression;

namespace SequencingFiles
{
    public abstract class GzipCommon : ClosableDisposable
    {
        #region member variables

        protected const byte LineFeedChar = 10;
        protected int CurrentOffset;
        protected long BufferOffset;
        protected string FilePath;
        protected IntPtr FileStreamPointer;
        protected bool IsOpen;
        protected byte[] LineBuffer;

        #endregion

        /// <summary>
        ///     constructor
        /// </summary>
        protected GzipCommon()
        {
            IsOpen = false;

            CurrentOffset = 0;
            BufferOffset = 0;
            FileStreamPointer = IntPtr.Zero;
        }

        /// <summary>
        /// returns the current position
        /// N.B. only tested with GzipReader so far
        /// </summary>
        public long GetCurrentPosition()
        {
            return BufferOffset + CurrentOffset;
        }

        /// <summary>
        ///     Opens the specified filename
        /// </summary>
        protected void Open(string filename, string fileMode)
        {
            // sanity check
            if (IsOpen)
            {
                throw new ApplicationException("ERROR: The Open method was called even though the file is already open.");
            }

            FileStreamPointer = SafeNativeMethods.gzopen(filename, fileMode);

            if (FileStreamPointer == IntPtr.Zero)
            {
                //for some reason throwing an exception here can hang. Possibly from the managed to native transition?
                //if you encounter that you should check for the complicating condition before trying to open the gzip file and handle it accordingly
                //you may have to call Environment.Exit(-1) as throwing an exception may hang
                //interestingly enabling native code debugging is a workaround but that won't work outside Visual Studio
                Console.Error.WriteLine(string.Format("ERROR: Unable to open the file ({0}) for reading", filename));
                Console.Error.WriteLine("This process will now exit");
                System.Diagnostics.StackTrace t = new System.Diagnostics.StackTrace();
                Console.Error.WriteLine(t.ToString());
                System.Environment.Exit(-1);
            }

            FilePath = filename;
            IsOpen = true;
            CurrentOffset = 0;
            BufferOffset = 0;
        }
    }
}