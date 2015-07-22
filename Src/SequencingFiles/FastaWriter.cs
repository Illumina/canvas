using System;
using System.IO;

namespace SequencingFiles
{
    public class FastaWriter : IDisposable
    {
        #region member variables

        private const int FastaLineLength = 70;
        private bool _mIsDisposed;
        private bool _mIsOpen;
        private StreamWriter _mWriter;
        public bool SkipLineWrapping { get; set; }

        #endregion

        // constructor
        public FastaWriter(string filename)
        {
            _mIsOpen = false;
            _mIsDisposed = false;
            SkipLineWrapping = false;
            Open(filename);
        }

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        // destructor
        ~FastaWriter()
        {
            Dispose(false);
        }

        /// <summary>
        ///     Closes the file
        /// </summary>
        public void Close()
        {
            if (_mIsOpen)
            {
                _mIsOpen = false;
                _mWriter.Close();
                _mWriter.Dispose();
            }
        }

        protected virtual void Dispose(bool disposing)
        {
            lock (this)
            {
                if (!_mIsDisposed)
                {
                    _mIsDisposed = true;
                    Close();
                }
            }
        }

        // Implement IDisposable

        /// <summary>
        ///     Opens the file
        /// </summary>
        public void Open(string filename)
        {
            _mIsOpen = true;
            _mWriter = new StreamWriter(filename) {NewLine = "\n"};
        }

        /// <summary>
        ///     Writes the FASTA entry
        /// </summary>
        public void WriteEntry(string sequenceName, string sequenceBases)
        {
            // write the header
            _mWriter.WriteLine(">{0}", sequenceName);

            // write the bases
            if (SkipLineWrapping)
            {
                _mWriter.WriteLine(sequenceBases);
            }
            else
            {
                char[] bases = sequenceBases.ToCharArray();
                int remainingBases = bases.Length;
                int offset = 0;

                while (remainingBases > FastaLineLength)
                {
                    _mWriter.WriteLine("{0}", new string(bases, offset, FastaLineLength));
                    offset += FastaLineLength;
                    remainingBases -= FastaLineLength;
                }

                if (remainingBases > 0)
                {
                    _mWriter.WriteLine("{0}", new string(bases, offset, remainingBases));
                }
            }
        }
    }
}