using System;
using System.IO;

namespace SequencingFiles
{
    public class IndexedFastaLoader
    {
        private readonly StreamReader _indexReader;
        private readonly StreamReader _fastaReader;

        public IndexedFastaLoader(StreamReader indexReader, StreamReader fastaReader)
        {
            _indexReader = indexReader;
            _fastaReader = fastaReader;
        }

        public string LoadFastaSequence(string chromosomeName)
        {
            AdvanceStreamToSequenceDescription(chromosomeName);
            return new FastaLoader(_fastaReader).LoadFastaSequence(chromosomeName);
        }

        private void AdvanceStreamToSequenceDescription(string chromosomeName)
        {
            long byteOffset = GetByteOffsetForSequence(chromosomeName);
            _fastaReader.TrueSeek(byteOffset, SeekOrigin.Begin);
        }

        private long GetByteOffsetForSequence(string chromosomeName)
        {
            _indexReader.TrueSeek(0, SeekOrigin.Begin);
            long byteOffset = 0;
            while (true)
            {
                string fileLine = _indexReader.ReadLine();
                if (fileLine == null) break;
                string[] bits = fileLine.Split('\t');
                if (bits[0] == chromosomeName)
                {
                    return byteOffset;
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
            throw new ArgumentException($"{chromosomeName} was not found in the FASTA index");
        }
    }
}