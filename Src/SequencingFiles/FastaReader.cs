using System;
using System.IO;
using System.Text;
using System.Text.RegularExpressions;

namespace SequencingFiles
{
    public class FastaReader : IDisposable
    {
        #region member variables
        private readonly Regex _mNameRegex;
        private bool _mIsDisposed;
        private bool _mIsOpen;
        private StreamReader _mReader;
        public bool SkipHeaderParsing { get; set; }
        #endregion

        /// <summary>
        ///     constructor
        /// </summary>
        /// <param name="filename"></param>
        public FastaReader(string filename)
        {
            _mNameRegex = new Regex(@"^>(\S+)", RegexOptions.Compiled);
            _mIsOpen = false;
            _mIsDisposed = false;
            SkipHeaderParsing = false;
            Open(filename);
        }

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        // destructor
        ~FastaReader()
        {
            Dispose(false);
        }

        public void Seek(long position)
        {
            _mReader.BaseStream.Seek(position, SeekOrigin.Begin);
        }

        /// <summary>
        ///     Closes the file
        /// </summary>
        public void Close()
        {
            if (_mIsOpen)
            {
                _mIsOpen = false;
                _mReader.Close();
                _mReader.Dispose();
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
        ///     Reads the FASTA entry
        /// </summary>
        /// <returns>Returns false if no more entries are available.</returns>
        public bool GetNextEntry(ref GenericRead sequence)
        {
            // read the header
            string header = _mReader.ReadLine();
            if (header == null) return false;

            // sanity check
            if (!header.StartsWith(">"))
            {
                throw new ApplicationException("Encountered a FASTA header that did not start with '>'");
            }

            // extract the sequence name
            if (SkipHeaderParsing)
            {
                sequence.Name = header.Substring(1);
            }
            else
            {
                Match nameMatch = _mNameRegex.Match(header);
                sequence.Name = nameMatch.Groups[1].Value;
            }

            // read the bases
            StringBuilder sb = new StringBuilder();
            int peek = _mReader.Peek();
            while ((peek != -1) && (peek != '>'))
            {
                string line = _mReader.ReadLine();
                if (line == null) break;
                sb.Append(line);
                peek = _mReader.Peek();
            }

            sequence.Bases = sb.ToString();

            return true;
        }

        /// <summary>
        ///     Opens the file
        /// </summary>
        public void Open(string filename)
        {
            _mIsOpen = true;
            _mReader = new StreamReader(filename);
        }
    }
}