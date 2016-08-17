using System;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using SequencingFiles;

namespace CanvasCommon
{
    public class CanvasIO
    {
        public static void WriteToTextFile(string outfile, IEnumerable<GenomicBin> bins)
        {
            using (GzipWriter writer = new GzipWriter(outfile))
            {
                foreach (GenomicBin bin in bins)
                {
                    writer.WriteLine(string.Format("{0}\t{1}\t{2}\t{3:F2}\t{4}", bin.Chromosome, bin.Start, bin.Stop, bin.Count, bin.GC));
                }
            }
        }

        public static IEnumerable<GenomicBin> IterateThroughTextFile(string infile)
        {
            using (GzipReader reader = new GzipReader(infile))
            {
                string row;

                while ((row = reader.ReadLine()) != null)
                {

                    string[] fields = row.Split('\t');

                    string chr = fields[0];
                    int start = Convert.ToInt32(fields[1]);
                    int stop = Convert.ToInt32(fields[2]);
                    float count = float.Parse(fields[3]);
                    int gc = Convert.ToInt32(fields[4]);

                    GenomicBin bin = new GenomicBin(chr, start, stop, gc, count);
                    yield return bin;
                }
            }
        }

        public static List<GenomicBin> ReadFromTextFile(string infile)
        {
            return IterateThroughTextFile(infile).ToList();
        }

        /// <summary>
        /// Gets GenomicBins by chromosome as an OrderedDictionary. The order of chromosomes in path is preserved.
        /// </summary>
        /// <param name="path"></param>
        /// <returns>OrderedDictionary with string key and List<GenomicBin> value</returns>
        public static OrderedDictionary GetGenomicBinsByChrom(string path)
        {
            OrderedDictionary binsByChrom = new OrderedDictionary();

            foreach (var bin in IterateThroughTextFile(path))
            {
                if (!binsByChrom.Contains(bin.Chromosome))
                {
                    binsByChrom[bin.Chromosome] = new List<GenomicBin>();
                }
                ((List<GenomicBin>)binsByChrom[bin.Chromosome]).Add(bin);
            }

            return binsByChrom;
        }

        // write localSD metric
        public static void WriteLocalSDToTextFile(string outfile, double localSD) 
        {
            using (StreamWriter writer = new StreamWriter(outfile))
            {
                writer.Write("#localSD\t" + localSD);
                writer.WriteLine();
            }       
        }

        // read localSD metric
        public static double ReadLocalSDFromTextFile(string infile)
        {
            double localSDmetric = -1.0;
            using (StreamReader reader = new StreamReader(infile))
            {
                string row;

                while ((row = reader.ReadLine()) != null)
                {
                    int tabs = (int)row.Count(ch => ch == '\t');
                    if (tabs > 0) 
                    {
                        string[] fields = row.Split('\t');
                        string localSDstring = fields[0];
                        if (localSDstring == "#localSD")
                            localSDmetric = Convert.ToDouble(fields[1]);                   
                    }
                }
            }

            return localSDmetric;
        }

        private static Dictionary<string, string> GetChromosomeAlternativeNames(IEnumerable<string> keys)
        {
            Dictionary<string, string> results = new Dictionary<string,string>();
            foreach (string key in keys)
            {
                if (key.StartsWith("chr"))
                {
                    results[key.Replace("chr", "")] = key;
                }
                else
                {
                    results["chr" + key] = key;
                }

            }
            return results;

        }

        /// <summary>
        /// Parse the outputs of CanvasSNV, and note these variant frequencies in the appropriate segment.
        /// </summary>
        public static float LoadVariantFrequencies(string variantFrequencyFile, List<CanvasSegment> segments)
        {
            
            Console.WriteLine("{0} Load variant frequencies from {1}", DateTime.Now, variantFrequencyFile);
            int count = 0;
            Dictionary<string, List<CanvasSegment>> segmentsByChromosome = CanvasSegment.GetSegmentsByChromosome(segments);
            Dictionary<string, string> alternativeNames = GetChromosomeAlternativeNames(segmentsByChromosome.Keys);
            long totalCoverage = 0;
            int totalRecords = 0;
            using (GzipReader reader = new GzipReader(variantFrequencyFile))
            {
                while (true)
                {
                    string fileLine = reader.ReadLine();
                    if (fileLine == null) break;
                    if (fileLine.Length == 0 || fileLine[0] == '#') continue; // Skip headers
                    string[] bits = fileLine.Split('\t');
                    if (bits.Length < 6)
                    {
                        Console.Error.WriteLine("* Bad line in {0}: '{1}'", variantFrequencyFile, fileLine);
                        continue;
                    }
                    string chromosome = bits[0];
                    if (!segmentsByChromosome.ContainsKey(chromosome))
                    {
                        if (alternativeNames.ContainsKey(chromosome))
                        {
                            chromosome = alternativeNames[chromosome];
                        }
                        else continue;
                    }

                    int position = int.Parse(bits[1]); // 1-based (from the input VCF to Canvas SNV)
                    int countRef = int.Parse(bits[4]);
                    int countAlt = int.Parse(bits[5]);
                    if (countRef + countAlt < 10) continue;
                    float VF = countAlt / (float)(countRef + countAlt);
                    // Binary search for the segment this variant hits:
                    List<CanvasSegment> chrSegments = segmentsByChromosome[chromosome];
                    int start = 0;
                    int end = chrSegments.Count - 1;
                    int mid = (start + end) / 2;
                    while (start <= end)
                    {
                        if (chrSegments[mid].End < position) // CanvasSegment.End is already 1-based
                        {
                            start = mid + 1;
                            mid = (start + end) / 2;
                            continue;
                        }
                        if (chrSegments[mid].Begin + 1 > position) // Convert CanvasSegment.Begin to 1-based by adding 1
                        {
                            end = mid - 1;
                            mid = (start + end) / 2;
                            continue;
                        }
                        chrSegments[mid].VariantFrequencies.Add(VF);
                        chrSegments[mid].VariantTotalCoverage.Add(countRef + countAlt);
                        count++;
                        totalCoverage += (countRef + countAlt); // use only coverage information in segments
                        totalRecords++;
                        break;
                    }
                }
            }
            float meanCoverage = 0;
            if (totalRecords > 0)
                meanCoverage = totalCoverage / Math.Max(1f, totalRecords);
            Console.WriteLine("{0} Loaded a total of {1} usable variant frequencies", DateTime.Now, count);
            return meanCoverage;
        }
    }
}
