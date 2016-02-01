using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace SequencingFiles
{
    public class VcfWriter : FileCommon
    {
        #region member variables
        public List<string> HeaderLines = new List<string>();
        public List<string> Samples = new List<string>();
        private BgzipWriter _writer;
        #endregion

        /// <summary>
        ///     constructor
        /// </summary>
        public VcfWriter(string filePath, List<string> headerLines = null, List<string> samples = null)
        {
            if (headerLines != null) { HeaderLines = headerLines; }
            if (samples != null) { Samples = samples; }
            Open(filePath);
        }

        public VcfWriter(string filePath, VcfReader reader) 
        {
            HeaderLines = reader.HeaderLines.ToList();
            Samples = reader.Samples.ToList();
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
            _writer = new BgzipWriter(filePath, GzipCompressionLevel.BestSpeed1);
            FileName = filePath;
            IsOpen = true;
            WriteHeaderLines();
        }

        private void WriteHeaderLines() 
        {
            foreach (string line in HeaderLines) 
            {
                if (line.StartsWith(VcfCommon.ChromosomeHeader)) { break; }
                _writer.WriteLine(line);
            }

            StringBuilder builder = new StringBuilder("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
            if (Samples.Any()) 
            {
                builder.Append("\tFORMAT");
                foreach (string sample in Samples) { builder.AppendFormat("\t{0}", sample); }
            }
            _writer.WriteLine(builder.ToString());
        }

        /// <summary>
        ///     Writes the variant
        /// </summary>
        public void WriteVariant(VcfVariant variant)
        {
            // sanity check
            if (!IsOpen)
            {
                throw new ApplicationException("ERROR: An attempt was made to write a variant to an unopened file.");
            }

            _writer.WriteLine(variant.ToString());
        }
    }
}