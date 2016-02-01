using System;
using System.Text.RegularExpressions;

namespace SequencingFiles
{
    public class FastqReader : FileCommon
    {
        #region member variables
        private GzipReader reader;
        private readonly Regex headerRegex;
        public bool SkipReadNameParsing { get; set; }
        //added as an option for reading SRA files with non standard headers
        public bool SkipAllReadNameParsing { get; set; }
        #endregion

        /// <summary>
        ///     constructor
        /// </summary>
        public FastqReader(string filename, bool skipReadNameParsing = false)
        {
            headerRegex = new Regex(@"^@(\S+)(?:\s+(\d):(\w):\d+:(\S+))?", RegexOptions.Compiled);
            SkipReadNameParsing = skipReadNameParsing;
            SkipAllReadNameParsing = false;
            Open(filename);
        }

        /// <summary>
        ///     Closes the file
        /// </summary>
        public override void Close()
        {
            if (IsOpen)
            {
                IsOpen = false;
                reader.Close();
                reader.Dispose();
            }
        }

        /// <summary>
        ///     Opens the file
        /// </summary>
        public override void Open(string filename)
        {
            IsOpen = true;
            reader = new GzipReader(filename);
        }

        /// <summary>
        ///   Reads a FASTQ entry
        ///   Note that this assumes the ILMN format (4 lines per read). FASTQ files with multiline 
        ///     sequences will likely cause issues
        /// </summary>
        /// <returns>Returns false if no more entries are available.</returns>
        public bool GetNextFastqEntry(ref BoltRead read)
        {
            // grab the next entry
            bool foundProblem = false;
            string header;
            if ((header = reader.ReadLine()) == null) return false;
            if ((read.Bases = reader.ReadLine()) == null) foundProblem = true;
            if ((reader.ReadLine()) == null) foundProblem = true;
            if ((read.Qualities = reader.ReadLine()) == null) foundProblem = true;

            if (foundProblem)
            {
                throw new ApplicationException(
                    string.Format("ERROR: Unable to read the entire FASTQ entry in {0}. Is the file truncated?",
                                  FileName));
            }

            read.Header = header;
            if (!SkipAllReadNameParsing)
            {
                // parse the secondary information
                Match headerMatch = headerRegex.Match(header);
                if (!headerMatch.Success)
                {
                    throw new ApplicationException(string.Format("Unexpected FastQ header {0}", header));
                }

                read.UnparsedName = headerMatch.Groups[1].Value;

                if (headerMatch.Groups.Count > 2)
                {
                    read.ReadNum = int.Parse(headerMatch.Groups[2].Value);
                    read.IsFiltered = (headerMatch.Groups[3].Value[0] == 'Y');
                    read.Index = headerMatch.Groups[4].Value;
                }

                // parse the read name
                if (!SkipReadNameParsing)
                {
                    if (headerMatch.Groups.Count != 5 || !read.ParseReadName())
                    {
                        throw new ApplicationException(string.Format("Invalid field count in FastQ header {0}", header));
                    }
                }
            }
            return true;
        }
    }
}