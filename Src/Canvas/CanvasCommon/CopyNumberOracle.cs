using System;
using System.Collections.Generic;
using System.IO;
using Isas.SequencingFiles;

namespace CanvasCommon
{

    public class CNInterval
    {
        public int Start;
        public int End;
        public int CN;
        public double Heterogeneity;
    }

    /// <summary>
    /// This class represents a set of known copy number calls ("truth set"), useful for training models.
    /// </summary>
    public class CopyNumberOracle
    {
        #region Members
        public Dictionary<string, List<CNInterval>> KnownCN;
        #endregion

        public void LoadKnownCN(string oraclePath)
        {
            if (!File.Exists(oraclePath))
            {
                throw new ArgumentException(string.Format("* Error: Truth vcf not found at '{0}'", oraclePath));
            }

            if (oraclePath.EndsWith(".bed"))
            {
                LoadKnownCNBed(oraclePath);
                return;
            }
            LoadKnownCNVCF(oraclePath);
        }

        /// <summary>
        /// Load known CN data from a .bed file.  File lines have fields:
        /// chromosome, start, end, chromcountA, chromcountB
        /// So, copy number is the sum of the last 2 fields, major chromosome count is the max of the last 2 fields.
        /// </summary>
        /// <param name="oracleBedPath"></param>
        protected void LoadKnownCNBed(string oracleBedPath)
        {
            bool stripChr = false;
            int count = 0;
            this.KnownCN = new Dictionary<string, List<CNInterval>>();
            using (FileStream stream = new FileStream(oracleBedPath, FileMode.Open, FileAccess.Read))
            using (StreamReader reader = new StreamReader(stream))
            {
                while (true)
                {
                    string fileLine = reader.ReadLine();
                    if (fileLine == null) break;
                    if (fileLine.Length == 0 || fileLine[0] == '#') continue;
                    string[] bits = fileLine.Split('\t');
                    string chromosome = bits[0];
                    if (stripChr) chromosome = chromosome.Replace("chr", "");
                    if (!KnownCN.ContainsKey(chromosome)) KnownCN[chromosome] = new List<CNInterval>();
                    CNInterval interval = new CNInterval();
                    interval.Start = int.Parse(bits[1]);
                    interval.End = int.Parse(bits[2]);
                    interval.CN = int.Parse(bits[3]) + int.Parse(bits[4]);
                    if (bits.Length > 5)
                        interval.Heterogeneity = double.Parse(bits[5]);
                    else
                        interval.Heterogeneity = -1.0;
                    KnownCN[chromosome].Add(interval);
                    count++;
                }
            }
            Console.WriteLine(">>>Loaded {0} known-CN intervals", count);
        }

        public int GetKnownCNForSegment(CanvasSegment segment)
        {
            // Handle switched chromosome naming convention transparently:
            string chr = segment.Chr;
            if (!this.KnownCN.ContainsKey(segment.Chr))
            {
                chr = segment.Chr.Replace("chr", "");
                if (!this.KnownCN.ContainsKey(chr))
                {
                    chr = "chr" + segment.Chr;
                    if (!this.KnownCN.ContainsKey(chr)) return -1;
                }
            }
            int CN = -1;
            foreach (CNInterval interval in this.KnownCN[chr])
            {
                if (interval.End < segment.Begin) continue;
                if (interval.Start > segment.End) continue;
                int start = Math.Max(segment.Begin, interval.Start);
                int end = Math.Min(segment.End, interval.End);
                if ((end - start) * 2 >= (segment.Length))
                {
                    CN = interval.CN;
                    break;
                }
            }
            return CN;
        }

        public double GetKnownClonalityForSegment(CanvasSegment segment)
        {
            // Handle switched chromosome naming convention transparently:
            string chr = segment.Chr;
            if (!this.KnownCN.ContainsKey(segment.Chr))
            {
                chr = segment.Chr.Replace("chr", "");
                if (!this.KnownCN.ContainsKey(chr))
                {
                    chr = "chr" + segment.Chr;
                    if (!this.KnownCN.ContainsKey(chr)) return -1;
                }
            }
            double Clonality = -1;
            foreach (CNInterval interval in this.KnownCN[chr])
            {
                if (interval.End < segment.Begin) continue;
                if (interval.Start > segment.End) continue;
                int start = Math.Max(segment.Begin, interval.Start);
                int end = Math.Min(segment.End, interval.End);
                if ((end - start) * 2 >= (segment.Length))
                {
                    Clonality = interval.Heterogeneity;
                    break;
                }
            }
            return Clonality;
        }

        protected void LoadKnownCNVCF(string oracleVCFPath)
        {
            bool stripChr = false;

            // Load our "oracle" of known copy numbers:
            this.KnownCN = new Dictionary<string, List<CNInterval>>();
            int count = 0;
            using (GzipReader reader = new GzipReader(oracleVCFPath))
            {
                while (true)
                {
                    string fileLine = reader.ReadLine();
                    if (fileLine == null) break;
                    if (fileLine.Length == 0 || fileLine[0] == '#') continue;
                    string[] bits = fileLine.Split('\t');
                    if (bits.Length == 1 && bits[0].Trim().Length == 0) continue; // skip empty lines!
                    string chromosome = bits[0];
                    if (stripChr) chromosome = chromosome.Replace("chr", "");
                    if (!KnownCN.ContainsKey(chromosome)) KnownCN[chromosome] = new List<CNInterval>();
                    CNInterval interval = new CNInterval();
                    interval.Start = int.Parse(bits[1]);
                    interval.CN = -1;
                    string[] infoBits = bits[7].Split(';');
                    foreach (string subBit in infoBits)
                    {
                        if (subBit.StartsWith("CN="))
                        {
                            float tempCN = float.Parse(subBit.Substring(3));
                            if (subBit.EndsWith(".5"))
                            {
                                interval.CN = (int)Math.Round(tempCN + 0.1); // round X.5 up to X+1
                            }
                            else
                            {
                                interval.CN = (int)Math.Round(tempCN); // Round off
                            }
                        }
                        if (subBit.StartsWith("END="))
                        {
                            interval.End = int.Parse(subBit.Substring(4));
                        }
                    }
                    // Parse CN from Canvas output:
                    if (bits.Length > 8)
                    {
                        string[] subBits = bits[8].Split(':');
                        string[] subBits2 = bits[9].Split(':');
                        for (int subBitIndex = 0; subBitIndex < subBits.Length; subBitIndex++)
                        {
                            if (subBits[subBitIndex] == "CN")
                            {
                                interval.CN = int.Parse(subBits2[subBitIndex]);
                            }
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

    }
}
