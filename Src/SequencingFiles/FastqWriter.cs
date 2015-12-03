using System;

namespace SequencingFiles
{
    public class FastqWriter : FileCommon
    {
        #region member variables

        private GzipWriter _writer;
        public bool Append = true;

        #endregion

        /// <summary>
        ///     constructor
        /// </summary>
        public FastqWriter(string filePath, bool append = true)
        {
            Append = append;
            Open(filePath);
        }

        /// <summary>
        ///     Closes the file
        /// </summary>
        public override void Close()
        {
            if (IsOpen)
            {
                IsOpen = false;
                _writer.Close();
                _writer.Dispose();
            }
        }

        /// <summary>
        ///     Opens the file
        /// </summary>
        public override void Open(string filePath)
        {
            _writer = new GzipWriter(filePath, GzipCompressionLevel.BestSpeed1, Append);
            FileName = filePath;
            IsOpen = true;
        }

        /// <summary>
        ///     Writes the FASTQ entry
        /// </summary>
        public void WriteFastqEntry(BoltRead read)
        {
            // sanity check
            if (!IsOpen)
            {
                throw new ApplicationException("ERROR: An attempt was made to write a FASTQ entry to an unopened file.");
            }
            if (read.Header != null)
            {
                _writer.WriteLine(read.Header);
            }
            else if (read.UMI != null)
            {
                _writer.WriteLine(string.Format("@{0}:{1}:{2}:{3}:{4}:{5}:{6}:{7} {8}:{9}:0:{10}",
                                               read.InstrumentName,
                                               read.RunID,
                                               read.FlowcellID,
                                               read.Lane,
                                               read.Tile,
                                               read.X,
                                               read.Y,
                                               read.UMI, 
                                               read.ReadNum,
                                               (read.IsFiltered ? 'Y' : 'N'),
                                               read.Index));
            }
            else
            {
                _writer.WriteLine(string.Format("@{0}:{1}:{2}:{3}:{4}:{5}:{6} {7}:{8}:0:{9}",
                                               read.InstrumentName,
                                               read.RunID,
                                               read.FlowcellID,
                                               read.Lane,
                                               read.Tile,
                                               read.X,
                                               read.Y,
                                               read.ReadNum,
                                               (read.IsFiltered ? 'Y' : 'N'),
                                               read.Index));
            }
            _writer.WriteLine(read.Bases);
            _writer.WriteLine("+");
            _writer.WriteLine(read.Qualities);
        }
    }
}