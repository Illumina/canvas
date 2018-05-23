using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Isas.Framework.Logging;
using Isas.SequencingFiles;
using Isas.SequencingFiles.Vcf;

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

        private const string LocalSdMetricName = "localSD";
        public static void WriteLocalSdMetricToTextFile(string filePath, double localSd)
        {
            WriteCoverageMetricToTextFile(filePath, LocalSdMetricName, localSd);
        }

        private const string EvennessMetricName = "evenness";
        public static void WriteEvennessMetricToTextFile(string filePath, double evenness)
        {
            WriteCoverageMetricToTextFile(filePath, EvennessMetricName, evenness);
        }

        private static void WriteCoverageMetricToTextFile(string filePath, string metricName, double metricValue)
        {
            File.WriteAllLines(filePath, $"#{metricName}\t{metricValue}".Yield());
        }

        public static double ReadLocalSdMetricFromTextFile(string filePath)
        {
            return ReadCoverageMetricFromTextFile(filePath, LocalSdMetricName);
        }

        public static double ReadEvennessMetricFromTextFile(string filePath)
        {
            return ReadCoverageMetricFromTextFile(filePath, EvennessMetricName);
        }

        private static double ReadCoverageMetricFromTextFile(string filePath, string metricName)
        {
            var lines = File.ReadAllLines(filePath);
            var metricLine = lines.Where(line => line.StartsWith($"#{metricName}")).ToList();
            if (!metricLine.Any())
                throw new ArgumentException($"Did not find {metricName} metric in file '{filePath}'");
            if (metricLine.Count > 1)
                throw new ArgumentException($"Found multiple {metricName} metrics in file '{filePath}'");

            var metricValue = metricLine.Single().Split('\t')[1];
            return double.Parse(metricValue);
        }

        public static Dictionary<string, List<Balleles>> ReadFrequenciesWrapper(ILogger logger,
            IFileLocation variantFrequencyFile, IReadOnlyDictionary<string, List<BedInterval>> intervalsByChromosome)
        {

            using (var reader = new GzipOrTextReader(variantFrequencyFile.FullName))
            {
                logger.Info($"Load variant frequencies from {variantFrequencyFile}");
                return ReadFrequencies(reader, intervalsByChromosome);
            }
        }

        public static Dictionary<string, List<Balleles>> ReadFrequencies(GzipOrTextReader variantFrequencyFileReader,
            IReadOnlyDictionary<string, List<BedInterval>> intervalByChromosome)
        {
            const int minCounts = 10;
            var alleleCountsByChromosome = new Dictionary<string, List<Balleles>>();
            foreach (string chr in intervalByChromosome.Keys)
                alleleCountsByChromosome[chr] = new List<Balleles>(intervalByChromosome[chr].Select(counter => new Balleles()));

            int index = 0;
            string prevChr = "";
            while (true)
            {
                string fileLine = variantFrequencyFileReader.ReadLine();
                if (fileLine == null) break;
                if (fileLine.Length == 0 || fileLine[0] == '#') continue; // Skip headers
                var columns = fileLine.Split('\t');
                string chr = columns[0];
                if (chr != prevChr)
                {
                    prevChr = chr;
                    index = 0;
                }
                int position = int.Parse(columns[1]); // 1-based (from the input VCF to Canvas SNV)
                int countRef = int.Parse(columns[4]);
                int countAlt = int.Parse(columns[5]);
                if (intervalByChromosome.Keys.All(chromosome => chromosome != chr))
                    continue;

                if (countRef + countAlt < minCounts) continue;

                // search forward for the right bin
                while (index < intervalByChromosome[chr].Count)
                {
                    if (intervalByChromosome[chr][index].End > position)
                        break;
                    ++index;
                }
                if (index >= intervalByChromosome[chr].Count)
                    continue;
                if (intervalByChromosome[chr][index].Start > position)
                    continue;

                alleleCountsByChromosome[chr][index].Add(new Ballele(position, countRef, countAlt));
            }
            return alleleCountsByChromosome;
        }
    }
}


