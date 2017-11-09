using System;
using System.Collections.Generic;
using System.IO;

namespace CanvasSNV
{
    /// <summary>
    /// Build up a histogram of variant frequency (VF) by region, or across all regions with a given copy number.
    /// Inputs:
    /// - Oracle VCF - gives the true copy number (wherever it's known) for the tumor genome
    /// - Folder of empirical variant frequencies, measured by CanvasSNV
    /// - Output path, for writing our results
    /// </summary>
    class HistogramVF
    {
        #region Members
        int[][] HistogramByCN = new int[10][];
        Dictionary<string, List<CNInterval>> KnownCN = new Dictionary<string, List<CNInterval>>();
        #endregion

        protected void LoadKnownCN(string oracleVCFPath)
        {
            // Load our "oracle" of known copy numbers:
            this.KnownCN = new Dictionary<string, List<CNInterval>>();
            int count = 0;
            using (FileStream stream = new FileStream(oracleVCFPath, FileMode.Open, FileAccess.Read))
            using (StreamReader reader = new StreamReader(stream))
            {
                while (true)
                {
                    string fileLine = reader.ReadLine();
                    if (fileLine == null) break;
                    if (fileLine.Length == 0 || fileLine[0] == '#') continue;
                    string[] bits = fileLine.Split('\t');
                    string chromosome = bits[0];
                    if (!KnownCN.ContainsKey(chromosome)) KnownCN[chromosome] = new List<CNInterval>();
                    CNInterval interval = new CNInterval();
                    interval.Start = int.Parse(bits[1]);
                    interval.CN = -1;
                    bits = bits[7].Split(';');
                    foreach (string subBit in bits)
                    {
                        if (subBit.StartsWith("CN="))
                        {
                            interval.CN = int.Parse(subBit.Substring(3));
                        }
                        if (subBit.StartsWith("END="))
                        {
                            interval.End = int.Parse(subBit.Substring(4));
                        }

                    }
                    if (interval.End == 0 || interval.CN < 0)
                    {
                        Console.WriteLine("Error - bogus record!");
                        Console.WriteLine(fileLine);
                    }
                    else
                    {
                        KnownCN[chromosome].Add(interval);
                        count++;
                    }
                }
            }
            Console.WriteLine(">>>Loaded {0} known-CN intervals", count);
        }

        protected void PopulateHistogramByCN(string empiricalVariantFrequencyFolder)
        {
            // Initialize histogram:
            int[][] HistogramByCN = new int[10][];
            for (int CN = 0; CN < HistogramByCN.Length; CN++)
            {
                HistogramByCN[CN] = new int[101];
            }

            // Read in the results files for each file, and use them to accumulate a histogram of variant frequency 
            // for each CN:
            foreach (string filePath in Directory.GetFiles(empiricalVariantFrequencyFolder))
            {
                string fileName = Path.GetFileName(filePath);
                if (!fileName.EndsWith("results.txt")) continue;
                using (FileStream stream = new FileStream(filePath, FileMode.Open, FileAccess.Read))
                using (StreamReader reader = new StreamReader(stream))
                {
                    while (true)
                    {
                        string fileLine = reader.ReadLine();
                        if (fileLine == null) break;
                        if (fileLine.Length == 0 || fileLine[0] == '#') continue;
                        string[] bits = fileLine.Split('\t');
                        string chromosome = bits[0];
                        int position = int.Parse(bits[1]);
                        int countRef = int.Parse(bits[4]);
                        int countAlt = int.Parse(bits[5]);
                        if (countRef + countAlt < 10) continue;
                        float VF = countAlt / (float)(countRef + countAlt);
                        int bin = (int)Math.Round(100 * VF);
                        int CN = -1;
                        if (!KnownCN.ContainsKey(chromosome))
                        {
                            Console.WriteLine("Warning: CN not known for {0}", chromosome);
                            continue;
                        }
                        foreach (CNInterval interval in KnownCN[chromosome])
                        {
                            if (interval.Start <= position && interval.End >= position)
                            {
                                CN = interval.CN;
                                break;
                            }
                        }
                        if (CN >= 0 && CN < HistogramByCN.Length)
                        {
                            HistogramByCN[CN][bin]++;
                        }
                    }
                }
            }
        }

        public int SummarizeStatsByRegion(string oracleVCFPath, string empiricalVariantFrequencyFolder, string outputPath)
        {
            this.LoadKnownCN(oracleVCFPath);

            // Read in the results files for each file, and use them to accumulate stats on variant frequency 
            // for each *region*:
            foreach (string filePath in Directory.GetFiles(empiricalVariantFrequencyFolder))
            {
                string fileName = Path.GetFileName(filePath);
                if (!fileName.EndsWith("results.txt")) continue;
                using (FileStream stream = new FileStream(filePath, FileMode.Open, FileAccess.Read))
                using (StreamReader reader = new StreamReader(stream))
                {
                    while (true)
                    {
                        string fileLine = reader.ReadLine();
                        if (fileLine == null) break;
                        if (fileLine.Length == 0 || fileLine[0] == '#') continue;
                        string[] bits = fileLine.Split('\t');
                        string chromosome = bits[0];
                        int position = int.Parse(bits[1]);
                        int countRef = int.Parse(bits[4]);
                        int countAlt = int.Parse(bits[5]);
                        if (countRef + countAlt < 10) continue;
                        float VF = countAlt / (float)(countRef + countAlt);
                        if (!KnownCN.ContainsKey(chromosome))
                        {
                            Console.WriteLine("Warning: CN not known for {0}", chromosome);
                            continue;
                        }
                        foreach (CNInterval interval in KnownCN[chromosome])
                        {
                            if (interval.Start <= position && interval.End >= position)
                            {
                                interval.Frequencies.Add(VF);
                                break;
                            }
                        }
                    }
                }
            }
            using (FileStream stream = new FileStream(outputPath, FileMode.Create, FileAccess.Write))
            using (StreamWriter writer = new StreamWriter(stream))
            {
                // Report variant frequency stats for each region:
                foreach (string chromosome in KnownCN.Keys)
                {
                    List<CNInterval> intervals = KnownCN[chromosome];
                    foreach (CNInterval interval in intervals)
                    {
                        // Write one histogram out for each region large enough to have a goodly number of variants:
                        int[] Histogram = new int[101];
                        if (interval.Frequencies.Count < 10000) continue;
                        int total = 0;
                        foreach (float VF in interval.Frequencies)
                        {
                            int bin = (int)Math.Round(VF * 100);
                            Histogram[bin]++;
                            total++;
                        }
                        writer.WriteLine();
                        writer.WriteLine("#{0}\t{1}\t{2}\t{3}\t", chromosome, interval.Start, interval.End, interval.CN);
                        for (int bin = 0; bin < Histogram.Length; bin++)
                        {
                            writer.WriteLine("{0}\t{1}\t{2}", bin, Histogram[bin], 100 * Histogram[bin] / (float)total);
                        }
                    }
                }
            }
            return 0;
        }

        public int BuildHistogramByCN(string oracleVCFPath, string empiricalVariantFrequencyFolder, string outputPath)
        {
            Console.WriteLine(">>> HistogramVF.Main()");
            this.LoadKnownCN(oracleVCFPath);
            this.PopulateHistogramByCN(empiricalVariantFrequencyFolder);

            // Summarize results:
            Console.WriteLine(">>> Summarize results to {0}", outputPath);
            using (FileStream stream = new FileStream(outputPath, FileMode.Create, FileAccess.Write))
            using (StreamWriter writer = new StreamWriter(stream))
            {
                writer.NewLine = "\n";
                writer.WriteLine("#Bin\tCN0\tCN1\tCN2\tCN3\tCN4\tCN5\tCN6\tCN7\tCN8\tCN9\t");
                for (int bin = 0; bin <= 100; bin++)
                {
                    writer.Write("{0}\t", bin);
                    for (int CN = 0; CN < HistogramByCN.Length; CN++)
                    {
                        writer.Write("{0}\t", HistogramByCN[CN][bin]);
                    }
                    writer.WriteLine();
                }
            }
            Console.WriteLine(">>> Results written to {0}", outputPath);

            return 0;
        }

        class CNInterval
        {
            public int Start;
            public int End;
            public int CN;
            public List<float> Frequencies = new List<float>();
        }
    }
}
