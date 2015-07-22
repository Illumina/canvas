using System;
using System.Collections.Generic;
using System.IO;
using System.Text;

namespace SequencingFiles
{
    public class DbSnpIndex
    {
        #region members

        private readonly Dictionary<string, long> _offsetIndex;

        #endregion

        // constructor
        public DbSnpIndex(string indexPath)
        {
            _offsetIndex = new Dictionary<string, long>();
            LoadIndex(indexPath);
        }

        /// <summary>
        ///     creates an index from a dbSNP file
        /// </summary>
        public static void CreateIndex(string dbSnpPath, string indexPath)
        {
            using (StreamReaderWithPosition reader = new StreamReaderWithPosition(dbSnpPath))
            using (StreamWriter writer = new StreamWriter(indexPath))
            {
                string oldRefSeqName = null;
                while (true)
                {
                    long filePosition = reader.Position;
                    string line = reader.ReadLine();
                    if (line == null) break;
                    if (line.Length == 0) continue;

                    // get the reference sequence name
                    string currentRefSeqName = GetReferenceSequenceName(line);

                    if (currentRefSeqName != oldRefSeqName)
                    {
                        writer.WriteLine("{0}\t{1}", currentRefSeqName, filePosition);
                        oldRefSeqName = currentRefSeqName;
                        Console.WriteLine("{0}\t{1}", currentRefSeqName, filePosition);
                    }
                }
            }
        }

        /// <summary>
        ///     Returns the reference sequence name given a line from the dbSNP file
        /// </summary>
        public static string GetReferenceSequenceName(string line)
        {
            if (string.IsNullOrEmpty(line)) return string.Empty;
            string[] columns = line.Split('\t');
            if (columns.Length < 2) return string.Empty;
            return columns[1];
        }

        /// <summary>
        ///     returns true if the file offset was found in the index, false otherwise
        /// </summary>
        public bool GetOffset(string chromosomeName, out long offset)
        {
            if (_offsetIndex.TryGetValue(chromosomeName, out offset)) return true;
            return false;
        }

        /// <summary>
        ///     loads the index from the specified file
        /// </summary>
        private void LoadIndex(string indexPath)
        {
            using (StreamReader reader = new StreamReader(indexPath))
            {
                while (true)
                {
                    // retrieve the next line
                    string line = reader.ReadLine();
                    if (string.IsNullOrEmpty(line)) break;

                    // store the information in our dictionary
                    string[] columns = line.Split('\t');

                    if (columns.Length != 2)
                    {
                        throw new ApplicationException(
                            string.Format("ERROR: Expected two columns in the dbSNP index, but found {0} columns",
                                          columns.Length));
                    }

                    long offset;
                    try
                    {
                        offset = long.Parse(columns[1]);
                    }
                    catch (Exception)
                    {
                        throw new ApplicationException(
                            string.Format(
                                "ERROR: Unable to convert the offset ({0} in the dbSNP index to a long integer.",
                                columns[1]));
                    }

                    _offsetIndex[columns[0]] = offset;
                }
            }
        }

        /// <summary>
        ///     Wrapper for a binary reader which enables reading an ASCII text file by line, and also
        ///     provides the underlying file byte-position.  (Note that streamreader.position is NOT a byte
        ///     position suitable for seeking to!)
        /// </summary>
        private class StreamReaderWithPosition : IDisposable
        {
            #region Members

            private readonly StringBuilder _lineBuilder;
            private byte[] _buffer;
            private int _bufferPos;
            private BinaryReader _inputReader;
            private Stream _inputStream;
            public long Position;

            #endregion

            public StreamReaderWithPosition(string filePath)
            {
                _bufferPos = 0;
                _buffer = new byte[0];
                //Buffer = new char[1024*1024];
                _inputStream = new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.Read);
                _inputReader = new BinaryReader(_inputStream);
                _lineBuilder = new StringBuilder();
                Position = 0;
            }

            public virtual void Dispose()
            {
                if (_inputReader != null)
                {
                    try
                    {
                        _inputReader.Close();
                    }
                    catch
                    {
                    }
                    _inputReader = null;
                }
                if (_inputStream != null)
                {
                    try
                    {
                        _inputStream.Close();
                    }
                    catch
                    {
                    }
                    _inputStream = null;
                }
            }

            public string ReadLine()
            {
                _lineBuilder.Length = 0;
                while (true)
                {
                    if (_bufferPos >= _buffer.Length)
                    {
                        // Read more:
                        _buffer = _inputReader.ReadBytes(1024*1024);
                        if (_buffer.Length == 0)
                        {
                            if (_lineBuilder.Length == 0) return null; // EOF!
                            return _lineBuilder.ToString();
                        }
                        _bufferPos = 0;
                    }
                    char letter = (char) _buffer[_bufferPos];
                    _bufferPos++;
                    Position++;
                    if (letter == '\n' || letter == '\r')
                    {
                        return _lineBuilder.ToString();
                    }
                    _lineBuilder.Append(letter);
                }
            }
        }
    }
}