using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace SequencingFiles
{
    public class BedEntry : IComparable<BedEntry>, IEquatable<BedEntry>
    {
        public string Chromosome { get; private set; }
        public int Start { get; private set; }
        public int End { get; private set; }
        public string Name { get; private set; }
        public int Score { get; private set; }
        public char Strand { get; private set; } // Defines the strand - either '+' or '-'.
        public int ThickStart { get; private set; }
        public int ThickEnd { get; private set; }
        public int ItemRgb { get; private set; }
        public int BlockCount { get; private set; }
        public ICollection<int> BlockSizes { get; private set; }
        public ICollection<int> BlockStarts { get; private set; }
        public override int GetHashCode()
        {
            return ToString().GetHashCode();
        }
        public override string ToString()
        {
            return string.Format("{0}:{1}:{2}", Chromosome, Start, End);
        }
        public override bool Equals(object obj)
        {
            var other = obj as BedEntry;
            if (other != null)
            {
                return Equals(other);
            }
            return base.Equals(obj);
        }

        public BedEntry(string line)
        {
            string[] parts = line.Split('\t');
            int column = 0;
            if (parts.Length < 3)
            {
                throw new Exception(string.Format("Invalid bed line '{0}' - Each line in the body of a .bed file must have at least 3 tab-delimited fields", line));
            }
            Chromosome = parts[column++];
            Start = Int32.Parse(parts[column++]);
            End = Int32.Parse(parts[column++]);
            if (column < parts.Length)
                Name = parts[column++];
            if (column < parts.Length)
                Score = Int32.Parse(parts[column++]);
            if (column < parts.Length)
                Strand = parts[column++][0] == '+' ? '+' : '-';
            if (column < parts.Length)
                ThickStart = Int32.Parse(parts[column++]);
            if (column < parts.Length)
                ThickEnd = Int32.Parse(parts[column++]);
            if (column < parts.Length)
                ItemRgb = Int32.Parse(parts[column++]);
            if (column < parts.Length)
                BlockCount = Int32.Parse(parts[column++]);
            if (column < parts.Length && BlockCount > 0)
                BlockSizes = parts[column++].Split(new[] { ',' }, StringSplitOptions.RemoveEmptyEntries).Select(x => Int32.Parse(x)).ToArray();
            if (column < parts.Length && BlockCount > 0)
                BlockStarts = parts[column++].Split(new[] { ',' }, StringSplitOptions.RemoveEmptyEntries).Select(x => Int32.Parse(x)).ToArray();
            if (BlockCount > 0 && (BlockSizes.Count != BlockCount || BlockStarts.Count != BlockCount))
                throw new InvalidDataException("Mismatch of BlockCount and elements");
        }

        public int CompareTo(BedEntry other)
        {
            if (object.ReferenceEquals(this, other))
                return 0;
            int result = Chromosome.CompareTo(other.Chromosome);
            if (0 == result)
                result = Start.CompareTo(other.Start);
            if (0 == result)
                result = End.CompareTo(other.End);
            return result;
        }

        public bool Equals(BedEntry other)
        {
            if (object.ReferenceEquals(this, other))
                return true;

            bool result = Chromosome == other.Chromosome && Start == other.Start && End == other.End &&
                Name == other.Name &&
                Score == other.Score &&
                Strand == other.Strand &&
                ThickStart == other.ThickStart &&
                ThickEnd == other.ThickEnd &&
                ItemRgb == other.ItemRgb &&
                BlockCount == other.BlockCount;
            if (result && BlockCount > 0)
                return BlockSizes.SequenceEqual(other.BlockSizes) && BlockStarts.SequenceEqual(other.BlockStarts);
            return result;
        }

        public ReferenceInterval CreateReferenceInterval()
        {
            return new ReferenceInterval(Chromosome, new Interval(Start + 1, End));
        }
    }
    public class BedReader : ClosableDisposable
    {
        #region Members
        TextReader _reader;
        public IEnumerable<string> HeaderLines { get; }
        private string BufferedLine;
        #endregion

        /// <summary>
        /// Parse headers specially - lines beginning with "#" or "track" or "browser" are headers:
        /// http://bedtools.readthedocs.org/en/latest/content/overview.html
        /// (Note: This .bed convention would create problems if you ever had a FASTA file whose first 
        /// contig started with #, track, or browser...)
        /// </summary>
        /// <param name="fileName"></param>
        public BedReader(string fileName) : this(File.OpenText(fileName))
        {
            var headerLines = new List<string>();
            while (true)
            {
                string fileLine = _reader.ReadLine();
                if (fileLine == null) break;
                if (fileLine.StartsWith("#") || fileLine.StartsWith("track") || fileLine.StartsWith("browser"))
                {
                    headerLines.Add(fileLine);
                    continue;
                }
                BufferedLine = fileLine;
                break;
            }
            HeaderLines = headerLines;
        }

        protected BedReader(TextReader input)
        {
            _reader = input;
        }

        public IEnumerable<BedEntry> GetEntries()
        {
            if (BufferedLine != null)
            {
                BedEntry entry = new BedEntry(BufferedLine);
                BufferedLine = null;
                yield return entry;
            }
            string line;
            while ((line = _reader.ReadLine()) != null)
                yield return new BedEntry(line);
        }

        public override void Close()
        {
            if (_reader != null)
            {
                _reader.Close();
                _reader = null;
            }
        }
    }
}
