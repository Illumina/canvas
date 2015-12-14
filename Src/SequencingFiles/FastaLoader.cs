using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace SequencingFiles
{
    public class FastaLoader
    {
        private readonly StreamReader _reader;
        private bool _nextSequenceHasEmptyDescription;

        public FastaLoader(StreamReader reader)
        {
            _reader = reader;
            InitializeStream();
        }

        private void InitializeStream()
        {
            // remove comment lines
            while (_reader.Peek() == ';') _reader.ReadLine();

            // the first sequence does not require a description line
            if (_reader.Peek() != '>')
                _nextSequenceHasEmptyDescription = true;
        }

        /// <summary>
        /// Not sure why anyone would want to do this, but this should maintain backwards compatibility
        /// </summary>
        public static string LoadFastaSequencesAndConcatenate(string fastaPath)
        {
            StringBuilder concatenatedSequence = new StringBuilder();
            foreach (var sequence in EnumerateFastaSequences(fastaPath))
            {
                concatenatedSequence.Append(sequence.Sequence);
            }
            return concatenatedSequence.ToString();
        }

        /// <summary>
        /// There can be only one Highlander...err Sequence
        /// </summary>
        public static string LoadSingleSequenceFasta(string fastaPath)
        {
            return EnumerateFastaSequences(fastaPath).Single().Sequence;
        }

        public static IEnumerable<FastaSequence> EnumerateFastaSequences(string fastaPath)
        {
            using (StreamReader reader = new StreamReader(fastaPath))
            {
                // make this an iterator so we don't prematurely dispose the reader
                foreach (var sequence in new FastaLoader(reader).EnumerateFastaSequences())
                    yield return sequence;
            }
        }

        public IEnumerable<FastaSequence> EnumerateFastaSequences()
        {
            while (true)
            {
                if (_reader.Peek() == -1 && !_nextSequenceHasEmptyDescription) yield break;
                yield return LoadNextFastaSequence(); ;
            }
        }

        public static Dictionary<string, string> LoadFastaSequences(string fastaPath)
        {
            using (StreamReader reader = new StreamReader(fastaPath))
            {
                return new FastaLoader(reader).LoadFastaSequences();
            }
        }

        public Dictionary<string, string> LoadFastaSequences()
        {
            var sequences = new Dictionary<string, string>();
            foreach (var sequence in EnumerateFastaSequences())
            {
                if (sequences.ContainsKey(sequence.Name))
                    throw new FileLoadException($"Found duplicated sequence {sequence.Name}");
                sequences[sequence.Name] = sequence.Sequence;
            }
            return sequences;
        }

        private FastaSequence LoadNextFastaSequence()
        {
            var chromosomeName = GetCurrentSequenceName();
            return new FastaSequence(chromosomeName, LoadCurrentSequence());
        }

        private string GetCurrentSequenceName()
        {
            if (_nextSequenceHasEmptyDescription)
            {
                _nextSequenceHasEmptyDescription = false;
                return "";
            }
            return GetSequenceFromLine(_reader.ReadLine());
        }

        public static string LoadFastaSequence(string fastaPath, string chromosomeName)
        {
            using (StreamReader fastaReader = new StreamReader(fastaPath))
            {
                string indexPath = string.Format("{0}.fai", fastaPath);
                if (!File.Exists(indexPath))
                    return new FastaLoader(fastaReader).LoadFastaSequence(chromosomeName);

                using (StreamReader indexReader = new StreamReader(indexPath))
                {
                    return new IndexedFastaLoader(indexReader, fastaReader).LoadFastaSequence(chromosomeName);
                }
            }
        }

        public string LoadFastaSequence(string chromosomeName)
        {
            AdvanceStreamToSequence(chromosomeName);
            return LoadCurrentSequence();
        }

        private string LoadCurrentSequence()
        {
            StringBuilder genomeBuilder = new StringBuilder();
            while (true)
            {
                int nextChar = _reader.Peek();
                if (nextChar == '>' || nextChar == -1) break; // leave the next sequence description on the stream
                genomeBuilder.Append(_reader.ReadLine().Trim());
            }
            return genomeBuilder.ToString();
        }

        private void AdvanceStreamToSequence(string chromosomeName)
        {
            if (TryAdvanceStreamToSequence(chromosomeName)) return;

            // try to search from the beginning of the stream
            _reader.TrueSeek(0, SeekOrigin.Begin);
            InitializeStream();
            if (!TryAdvanceStreamToSequence(chromosomeName))
                throw new ArgumentException($"Unable to find {chromosomeName} in fasta file");
        }

        private bool TryAdvanceStreamToSequence(string chromosomeName)
        {
            if (_nextSequenceHasEmptyDescription && chromosomeName == "") return true;
            while (true)
            {
                string line = _reader.ReadLine();
                if (line == null) break;
                if (line.FirstOrDefault() == '>' && GetSequenceFromLine(line) == chromosomeName)
                    return true;
            }
            return false;
        }

        private static string GetSequenceFromLine(string line)
        {
            if (line == null)
                throw new ApplicationException("FASTA file has no more sequences");
            return line.Substring(1).Split()[0];
        }
    }

    public class FastaSequence
    {
        public FastaSequence(string name, string sequence)
        {
            Name = name;
            Sequence = sequence;
        }

        public string Name { get; }
        public string Sequence { get; }
    }
}