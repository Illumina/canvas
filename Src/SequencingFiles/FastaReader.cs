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

        /// <summary>
        /// Load bases from a reference FASTA file.  If chromosomeName is set, load just that chromosome.
        /// If the .fai index is available, then use that to seek to the chromosome of interest, to save time!
        /// </summary>
        /// <returns></returns>
        public static string LoadFASTA(string fastaPath, string chromosomeName = null)
        {
            string faiPath = string.Format("{0}.fai", fastaPath);
            long? byteOffset = null;

            // Jump straight to the chromosome of interest:
            if (!string.IsNullOrEmpty(chromosomeName) && File.Exists(faiPath))
            {
                using (StreamReader reader = new StreamReader(faiPath))
                {
                    while (true)
                    {
                        string fileLine = reader.ReadLine();
                        if (fileLine == null) break;
                        string[] bits = fileLine.Split('\t');
                        if (bits[0] == chromosomeName)
                        {
                            break;
                        }
                        int lineCharacters = int.Parse(bits[3]);
                        int lineBytes = int.Parse(bits[4]);
                        long chrBases = long.Parse(bits[1]);
                        long thisChromosomeEnd = long.Parse(bits[2]);
                        thisChromosomeEnd += chrBases;
                        thisChromosomeEnd += (chrBases / lineCharacters) * (lineBytes - lineCharacters);
                        if (lineCharacters % chrBases != 0) thisChromosomeEnd += (lineBytes - lineCharacters);
                        byteOffset = thisChromosomeEnd;
                    }
                }
            }

            StringBuilder genomeBuilder = new StringBuilder();
            using (StreamReader reader = new StreamReader(fastaPath))
            {
                bool readThisChromosome = false;
                if (byteOffset != null) reader.TrueSeek((long)byteOffset, SeekOrigin.Begin);
                while (true)
                {
                    string fileLine = reader.ReadLine();
                    if (fileLine == null) break;
                    if (fileLine.Length == 0 || fileLine[0] == '>')
                    {
                        // We reached a chromosome name.
                        string name = fileLine.Substring(1).Split()[0];
                        readThisChromosome = false;
                        if (string.IsNullOrEmpty(chromosomeName))
                        {
                            if (genomeBuilder.Length > 0) throw new Exception(string.Format("Error: LoadFASTA called without a chromosome name, but FASTA file {0} has more than one contig", fastaPath));
                            readThisChromosome = true;
                        }
                        else if (string.Compare(name, chromosomeName, true) == 0) readThisChromosome = true;
                        if (genomeBuilder.Length > 0 && !readThisChromosome) break; // We already found our target.
                        continue;
                    }
                    if (readThisChromosome)
                        genomeBuilder.Append(fileLine.Trim());
                }
            }
            return genomeBuilder.ToString();
        }


    }
}