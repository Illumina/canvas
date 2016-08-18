using System;
using Isas.SequencingFiles;
using Isas.Shared.Utilities;

namespace FlagUniqueKmers
{
    /// <summary>
    /// Developer debug tool: Load the legacy 35-uniqueness genome.fa file, and a new output, and compare them.  
    /// How concordant are they?
    /// </summary>
    class CheckFlags
    {
        static public void CheckUniqueness()
        {
            string fastaPath = @"D:\Genomes\Homo_sapiens\UCSC\hg19\Sequence\WholeGenomeFasta\genome.fa";
            string[] Reads = new string[]
            {
                "AACCCTAACCCAACCCTAACCCTAACCCTAACCCT", // 10097 B
                "ACCCTAACCCAACCCTAACCCTAACCCTAACCCTA", // 10098 B
                "AGAGGACAACGCAGCTCCGCCCTCGCGGTGCTCTC", // 10553 A
                "TTTTTTCCTATACATACATACCCATGATAAAGTTT"  // 30763880 A
            };

            using (FastaReader readerA = new FastaReader(fastaPath))
            {
                GenericRead chrA = new GenericRead();
                while (true)
                {
                    bool result = readerA.GetNextEntry(ref chrA);
                    if (!result) break;
                    Console.WriteLine(chrA.Name);
                    string bases = chrA.Bases.ToUpperInvariant();
                    // Search for each:

                    for (int readIndex = 0; readIndex < Reads.Length; readIndex++)
                    {
                        int pos = -1;
                        while (true)
                        {
                            pos = bases.IndexOf(Reads[readIndex], pos + 1);
                            if (pos == -1) break;
                            Console.WriteLine("{0}\t{1}\t{2}\t{3}", readIndex, Reads[readIndex], chrA.Name, pos);
                        }
                        pos = -1;
                        string revComp = Utilities.GetReverseComplement(Reads[readIndex]);
                        while (true)
                        {
                            pos = bases.IndexOf(revComp, pos + 1);
                            if (pos == -1) break;
                            Console.WriteLine("{0}\t{1}\t{2}\t{3}\tRevComp", readIndex, Reads[readIndex], chrA.Name, pos);
                        }

                    }
                }
            }
            Console.WriteLine(">>>Done.");
        }

        static public int ProcessReferenceFASTA(string fastaPathA, string fastaPathB)
        {
            GenericRead chrA = new GenericRead();
            GenericRead chrB = new GenericRead();
            long CountAB = 0;
            long CountA = 0;
            long CountB = 0;
            long CountNeither = 0;
            using (FastaReader readerA = new FastaReader(fastaPathA))
            using (FastaReader readerB = new FastaReader(fastaPathB))
            {
                readerA.GetNextEntry(ref chrA); // Discard chrM from new output
                while (true)
                {
                    bool result = readerA.GetNextEntry(ref chrA);
                    if (!result) break;
                    readerB.GetNextEntry(ref chrB);
                    if (chrA.Bases.Length != chrB.Bases.Length) throw new Exception();
                    for (int baseIndex = 0; baseIndex < chrA.Bases.Length; baseIndex++)
                    {
                        bool isUniqueA = chrA.Bases[baseIndex] < 'a';
                        bool isUniqueB = chrB.Bases[baseIndex] < 'a';
                        if (isUniqueA && isUniqueB)
                        {
                            CountAB++;
                        }
                        else if (isUniqueA && !isUniqueB)
                        {
                            CountA++;
                        }
                        else if (!isUniqueA && isUniqueB)
                        {
                            CountB++;
                        }
                        else
                        {
                            CountNeither++;
                        }
                    }
                    Console.WriteLine("After {0}: {1},{2},{3},{4}", chrA.Name,
                        CountAB, CountA, CountB, CountNeither);
                    double percentAgreement = 100 * (CountAB + CountNeither) / (double)(CountAB + CountA + CountB + CountNeither);
                    Console.WriteLine("Percent agreement: {0:F2}", percentAgreement);
                }

            }
            return 0;
        }
    }
}
