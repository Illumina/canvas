using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using Illumina.Common;
using Isas.SequencingFiles;

namespace CanvasCommon
{
    public class CanvasIO
    {
        public static void WriteToTextFile(string outfile, IEnumerable<SampleGenomicBin> bins)
        {
            using (GzipWriter writer = new GzipWriter(outfile))
            {
                foreach (SampleGenomicBin bin in bins)
                {
                    writer.WriteLine(string.Format("{0}\t{1}\t{2}\t{3:F2}\t{4}", bin.GenomicBin.Chromosome, bin.Start, bin.Stop, bin.Count, bin.GenomicBin.GC));
                }
            }
        }

        public static IEnumerable<SampleGenomicBin> IterateThroughTextFile(string infile)
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

                    SampleGenomicBin bin = new SampleGenomicBin(chr, start, stop, gc, count);
                    yield return bin;
                }
            }
        }

        public static List<SampleGenomicBin> ReadFromTextFile(string infile)
        {
            return IterateThroughTextFile(infile).ToList();
        }

        /// <summary>
        /// Gets GenomicBins by chromosome as an OrderedDictionary. The order of chromosomes in path is preserved.
        /// Assumes that bins are sorted in ascending order by the start position within each chromosome.
        /// </summary>
        /// <param name="path"></param>
        /// <returns>OrderedDictionary with string key and List<GenomicBin> value</returns>
        public static OrderedDictionary<string, List<SampleGenomicBin>> GetGenomicBinsByChrom(string path)
        {
            OrderedDictionary<string, List<SampleGenomicBin>> binsByChrom = new OrderedDictionary<string, List<SampleGenomicBin>>();

            SampleGenomicBin prevBin = null;
            foreach (var bin in IterateThroughTextFile(path))
            {
                if (!binsByChrom.ContainsKey(bin.GenomicBin.Chromosome))
                {
                    binsByChrom[bin.GenomicBin.Chromosome] = new List<SampleGenomicBin>();
                    prevBin = null;
                }
                if (prevBin != null && bin.Start < prevBin.Start)
                    throw new Exception("Bins are not sorted in ascending order by the start position." + 
                        $" First offending bin: {bin.GenomicBin.Chromosome}\t{bin.Start}\t{bin.Stop}");

                binsByChrom[bin.GenomicBin.Chromosome].Add(bin);
                prevBin = bin;
            }

            return binsByChrom;
        }

        // write localSD metric
        public static void WriteLocalSDToTextFile(string outfile, double localSD) 
        {
            using (FileStream stream = new FileStream(outfile, FileMode.Create, FileAccess.Write))
            using (StreamWriter writer = new StreamWriter(stream))
            {
                writer.Write("#localSD\t" + localSD);
                writer.WriteLine();
            }       
        }

        // read localSD metric
        public static double ReadLocalSDFromTextFile(string infile)
        {
            double localSDmetric = -1.0;
            using (FileStream stream = new FileStream(infile, FileMode.Open, FileAccess.Read))
            using (StreamReader reader = new StreamReader(stream))
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

        public static Dictionary<string, string> GetChromosomeAlternativeNames(IEnumerable<string> keys)
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
        public static float LoadFrequenciesBySegment(string variantFrequencyFile, List<CanvasSegment> segments)
        {
            
            Console.WriteLine("{0} Load variant frequencies from {1}", DateTime.Now, variantFrequencyFile);
            int count = 0;
            var segmentsByChromosome = CanvasSegment.GetSegmentsByChromosome(segments);
            var alternativeNames = GetChromosomeAlternativeNames(segmentsByChromosome.Keys);
            int totalRecords;
            var intervalsByChromosome = new Dictionary<string, List<Interval>>();
            foreach (var chr in segmentsByChromosome.Keys)
            {
                intervalsByChromosome[chr] = new List<Interval>();
                foreach (var canvasSegment in segmentsByChromosome[chr])
                {
                    intervalsByChromosome[chr].Add(new Interval(canvasSegment.Begin, canvasSegment.End));
                }
            }
            var alleleCountsByChromosome = ReadFrequencies(variantFrequencyFile, intervalsByChromosome, alternativeNames, 
                out totalRecords, ref count);

            foreach (var chr in segmentsByChromosome.Keys)
            {
                for (int index = 0; index < segmentsByChromosome[chr].Count; index++)
                {
                    segmentsByChromosome[chr][index].Alleles.Frequencies.Add(countAlt / (float)(countRef + countAlt));
                    segmentsByChromosome[chr][index].Alleles.TotalCoverage.Add(countRef + countAlt);
                    segmentsByChromosome[chr][index].Alleles.Counts.Add(new Tuple<int, int>(countRef, countAlt));
                }
            }
            float meanCoverage = 0;
            if (totalRecords > 0)
                meanCoverage = totalCoverage / Math.Max(1f, totalRecords);
            Console.WriteLine("{0} Loaded a total of {1} usable variant frequencies", DateTime.Now, count);
            return meanCoverage;
        }

        private static long ReadFrequencies(string variantFrequencyFile, Dictionary<string, List<Interval>> intervalByChromosome,
            Dictionary<string, string> alternativeNames, out int totalRecords)
        {
            long totalCoverage = 0;
            int count;
            totalRecords = 0;
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
                    string chr = bits[0];
                    if (!intervalByChromosome.ContainsKey(chr))
                    {
                        if (alternativeNames.ContainsKey(chr))
                        {
                            chr = alternativeNames[chr];
                        }
                        else continue;
                    }

                    int position = int.Parse(bits[1]); // 1-based (from the input VCF to Canvas SNV)
                    int countRef = int.Parse(bits[4]);
                    int countAlt = int.Parse(bits[5]);
                    if (countRef + countAlt < 10) continue;
                    float VF = countAlt / (float) (countRef + countAlt);
                    // Binary search for the segment this variant hits:
                    int start = 0;
                    int end = intervalByChromosome[chr].Count - 1;
                    int mid = (start + end) / 2;
                    while (start <= end)
                    {
                        if (intervalByChromosome[chr][mid].OneBasedEnd < position) // CanvasSegment.End is already 1-based
                        {
                            start = mid + 1;
                            mid = (start + end) / 2;
                            continue;
                        }
                        if (intervalByChromosome[chr][mid].OneBasedStart + 1 > position) // Convert CanvasSegment.Begin to 1-based by adding 1
                        {
                            end = mid - 1;
                            mid = (start + end) / 2;
                            continue;
                        }

                        count++;
                        totalCoverage += countRef + countAlt; // use only coverage information in segments
                        totalRecords++;
                        break;
                    }
                }
            }
            return totalCoverage;
        }
    }
}
