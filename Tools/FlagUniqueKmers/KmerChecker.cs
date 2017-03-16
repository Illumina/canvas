using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using Isas.SequencingFiles;

namespace FlagUniqueKmers
{
    class KmerChecker
    {
        #region Members
        static private int KmerLength = 35;
        static private int MaxDictEntries = 400000000;  // We can go up this high, given the updates to app.config.
        Dictionary<string, long> Kmers = new Dictionary<string, long>(); // Dictionary: Kmer string -> first occurrence of the kmer in the genome
        List<BitArray> ChromosomeNonUniqueFlags = new List<BitArray>();
        List<BitArray> ChromosomeFinishedFlags = new List<BitArray>();
        private int PassIndex;
        private long GenomePosition;
        private long IncompletePositions;
        StringBuilder keyBuilder = new StringBuilder();
        #endregion

        /// <summary>
        /// Return the key for this kmer:
        /// - null, if it contains Ns
        /// - otherwise, compress both the string and its reverse complement to a short (more memory-efficient) string,
        ///   and return the first of the two.
        /// </summary>
        private string GetKeyForKmer(string kmer)
        {

            // Squish the kmer to a more byte-efficient representation:
            // 35 bases at 1 byte per character = 280 bits
            // using 2 bits per base, 35 bases can be encoded as 35*2=70 bits, which fits in 9 bytes.
            keyBuilder.Clear();
            byte currentChar = 0;
            bool badChar = false;
            for (int tempChar = 0; tempChar < kmer.Length; tempChar++)
            {
                currentChar *= 4;
                switch (kmer[tempChar])
                {
                    case 'A':
                        currentChar += 0;
                        break;
                    case 'C':
                        currentChar += 1;
                        break;
                    case 'G':
                        currentChar += 2;
                        break;
                    case 'T':
                        currentChar += 3;
                        break;
                    default:
                        badChar = true;
                        break;
                }
                if (tempChar % 4 == 3)
                {
                    keyBuilder.Append((char)currentChar);
                    currentChar = 0;
                }
            }
            keyBuilder.Append((char)currentChar);
            if (badChar) return null;
            string key = keyBuilder.ToString();

            // Now build key for the rev-comp:
            keyBuilder.Clear();
            currentChar = 0;
            int charCount = 0;
            for (int tempChar = kmer.Length - 1; tempChar >= 0; tempChar--, charCount++)
            {
                currentChar *= 4;
                switch (kmer[tempChar])
                {
                    case 'T':
                        currentChar += 0;
                        break;
                    case 'G':
                        currentChar += 1;
                        break;
                    case 'C':
                        currentChar += 2;
                        break;
                    case 'A':
                        currentChar += 3;
                        break;
                    default:
                        badChar = true;
                        break;
                }
                if (charCount % 4 == 3)
                {
                    keyBuilder.Append((char)currentChar);
                    currentChar = 0;
                }
            }
            keyBuilder.Append((char)currentChar);
            string key2 = keyBuilder.ToString();
            if (string.Compare(key, key2, StringComparison.Ordinal) < 0) return key;
            return key2;
        }

        
        private void ProcessOneChromosome(GenericRead fastaEntry, int chromosomeIndex)
        {
            BitArray nonUniqueFlags;
            BitArray finishedFlags;
            if (chromosomeIndex >= ChromosomeNonUniqueFlags.Count)
            {
                nonUniqueFlags = new BitArray(fastaEntry.Bases.Length);
                ChromosomeNonUniqueFlags.Add(nonUniqueFlags);
                finishedFlags = new BitArray(fastaEntry.Bases.Length);
                ChromosomeFinishedFlags.Add(finishedFlags);
            }
            nonUniqueFlags = ChromosomeNonUniqueFlags[chromosomeIndex];
            finishedFlags = ChromosomeFinishedFlags[chromosomeIndex];

            StringBuilder keyBuilder = new StringBuilder();

            string bases = fastaEntry.Bases.ToUpperInvariant();
            for (int startPos = 0; startPos < bases.Length; startPos++, GenomePosition++)
            {
                if (startPos % 1000000 == 0)
                {
                    Console.WriteLine(">>>{0} {1} {2} dict {3} incomplete {4}", PassIndex, fastaEntry.Name, startPos, Kmers.Keys.Count, IncompletePositions);
                }

                // Skip positions processed to completion in an earlier pass:
                if (finishedFlags[startPos]) continue;

                // Handle positions at very end of chromosome:
                if (startPos + KmerLength >= bases.Length)
                {
                    nonUniqueFlags.Set(startPos, true);
                    finishedFlags.Set(startPos, true);
                    continue;
                }

                // This position isn't completed yet.  Check its kmer against the dictionary:
                string kmer = bases.Substring(startPos, KmerLength);
                string key = GetKeyForKmer(kmer);

                if (key == null)
                {

                    // This position isn't a valid unique 35mer, because it has an N or some other bogus character.
                    nonUniqueFlags.Set(startPos, true);
                    finishedFlags.Set(startPos, true);
                    continue;
                }

                if (Kmers.ContainsKey(key))
                {

                    // This position is not unique.  Flag it as known non-unique:
                    nonUniqueFlags.Set(startPos, true);
                    finishedFlags.Set(startPos, true);

                    long oldPos = Kmers[key];
                    if (oldPos >= 0)
                    {
                        // Go back and flag the kmer position, even if on an old chromosome:
                        long tempPos = 0;
                        for (int tempIndex = 0; tempIndex < ChromosomeNonUniqueFlags.Count; tempIndex++)
                        {
                            if (tempPos + ChromosomeNonUniqueFlags[tempIndex].Length > oldPos)
                            {
                                // Note: Assumption is that OldPos - tempIndex will always be an int 
                                // since no single chromosome's length is longer than maxint.
                                int chrPos = (int)(oldPos - tempPos);
                                if (ChromosomeFinishedFlags[tempIndex][chrPos])
                                {
                                    throw new Exception("Error: Flagging an already-done position!");
                                }

                                ChromosomeNonUniqueFlags[tempIndex].Set(chrPos, true);
                                ChromosomeFinishedFlags[tempIndex].Set(chrPos, true);
                                break;
                            }
                            tempPos += ChromosomeNonUniqueFlags[tempIndex].Length;
                        }
                        Kmers[key] = -1;
                    }
                }
                else
                {
                    if (Kmers.Keys.Count >= MaxDictEntries)
                    {
                        IncompletePositions++;
                    }
                    else
                    {
                        Kmers[key] = GenomePosition;
                    }
                }
            } // loop over start positions
        }

        private void WriteOutputs(string fastaPath, string outputPath)
        {
            int chromosomeIndex = -1;
            using (FastaReader reader = new FastaReader(fastaPath))
            using (FastaWriter writer = new FastaWriter(outputPath))
            {
                GenericRead fastaEntry = new GenericRead();
                while (reader.GetNextEntry(ref fastaEntry))
                {
                    chromosomeIndex++;
                    StringBuilder baseBuilder = new StringBuilder();
                    BitArray nonUniqueFlags = ChromosomeNonUniqueFlags[chromosomeIndex];
                    for (int chromPos = 0; chromPos < fastaEntry.Bases.Length; chromPos++)
                    {
                        if (nonUniqueFlags[chromPos])
                        {
                            baseBuilder.Append(char.ToLowerInvariant(fastaEntry.Bases[chromPos]));
                        }
                        else
                        {
                            baseBuilder.Append(char.ToUpperInvariant(fastaEntry.Bases[chromPos]));
                        }
                    }
                    writer.WriteEntry(fastaEntry.Name, baseBuilder.ToString());
                }
            }
        }

        public void Main(string fastaPath, string outputPath)
        {
            Console.WriteLine("{0} Start", DateTime.Now);
            Console.WriteLine("Load FASTA file at {0}, write kmer-flagged output to {1}", fastaPath, outputPath);
            // Make multiple passes over the file, to manage memory usage.  In each pass, we'll accumulate a dictionary
            // of up to N kmers, and we'll flag as non-unique all occurrences of those kmers.  In subsequent passes,
            // we know we can ignore positions that were already flagged as non-unique.  We know that we're complete when
            // we make a pass where the dictionary doesn't fill up.
            this.PassIndex = 0;
            this.IncompletePositions = 1; // This will be set to the number of positions whose uniqueness status is UNKNOWN.
            HashSet<string> finishedChromosomes = new HashSet<string>(); // For speed, keep track of chromosomes that are already fully processed.
            while (IncompletePositions > 0)
            {
                IncompletePositions = 0;
                PassIndex++;
                Kmers.Clear();
                this.GenomePosition = 0;
                int chromosomeIndex = -1;
                using (FastaReader reader = new FastaReader(fastaPath))
                {
                    GenericRead fastaEntry = new GenericRead();
                    while (reader.GetNextEntry(ref fastaEntry))
                    {
                        chromosomeIndex++;
                        // Speedup option: Skip over chromosome that's already been fully processed:
                        if (finishedChromosomes.Contains(fastaEntry.Name))
                        {
                            GenomePosition += fastaEntry.Bases.Length;
                            continue;
                        }
                        ProcessOneChromosome(fastaEntry, chromosomeIndex);
                        if (IncompletePositions == 0 && !finishedChromosomes.Contains(fastaEntry.Name))
                        {
                            finishedChromosomes.Add(fastaEntry.Name);
                        }
                    } // loop over chromosomes

                    // Now let's go back and make a note of the various positions that are now definitively known
                    // to be unique kmers!
                    int UniqueCount = 0;
                    foreach (string key in Kmers.Keys)
                    {
                        long oldPos = Kmers[key];
                        if (oldPos < 0) continue;
                        UniqueCount++;
                        int tempPos = 0;
                        for (int tempIndex = 0; tempIndex < ChromosomeNonUniqueFlags.Count; tempIndex++)
                        {
                            if (tempPos + ChromosomeNonUniqueFlags[tempIndex].Length > oldPos)
                            {
                                // Note: Assumption is that OldPos - tempIndex will always be an int 
                                // since no single chromosome's length is longer than maxint.
                                ChromosomeFinishedFlags[tempIndex].Set((int)(oldPos - tempPos), true);
                                break;
                            }
                            tempPos += ChromosomeNonUniqueFlags[tempIndex].Length;
                        }
                    }
                    Console.WriteLine("{0} >>>Pass {1}: flagged {2} unique kmers", DateTime.Now, PassIndex, UniqueCount);
                    Console.WriteLine();

                }
            } // Pass loop

            Console.WriteLine("{0} Flagging complete", DateTime.Now);
            this.WriteOutputs(fastaPath, outputPath);
            Console.WriteLine("{0} Output written to {1}", DateTime.Now, outputPath);
        }
    }
}
