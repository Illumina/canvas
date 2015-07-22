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
            var parts = line.Split('\t');
            int i = 0;
            Chromosome = parts[i++];
            Start = Int32.Parse(parts[i++]);
            End = Int32.Parse(parts[i++]);
            if (i < parts.Length)
                Name = parts[i++];
            if (i < parts.Length)
                Score = Int32.Parse(parts[i++]);
            if (i < parts.Length)
                Strand = parts[i++][0] == '+' ? '+' : '-';
            if (i < parts.Length)
                ThickStart = Int32.Parse(parts[i++]);
            if (i < parts.Length)
                ThickEnd = Int32.Parse(parts[i++]);
            if (i < parts.Length)
                ItemRgb = Int32.Parse(parts[i++]);
            if (i < parts.Length)
                BlockCount = Int32.Parse(parts[i++]);
            if (i < parts.Length && BlockCount > 0)
                BlockSizes = parts[i++].Split(new [] {','}, StringSplitOptions.RemoveEmptyEntries).Select(x => Int32.Parse(x)).ToArray();
            if (i < parts.Length && BlockCount > 0)
                BlockStarts = parts[i++].Split(new[] { ',' }, StringSplitOptions.RemoveEmptyEntries).Select(x => Int32.Parse(x)).ToArray();
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
    }
    public class BedReader : ClosableDisposable
    {
        TextReader _reader;
        public BedReader(string fileName) : this(File.OpenText(fileName))
        { }
        public BedReader(TextReader input)
        {
            _reader = input;
        }
        public IEnumerable<BedEntry> GetEntries()
        {
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
