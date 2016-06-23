using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using SequencingFiles;

namespace CanvasCommon
{
    public class PloidyInfo
    {
        #region Members
        public string HeaderLine;
        public Dictionary<string, List<PloidyInterval>> PloidyByChromosome = new Dictionary<string, List<PloidyInterval>>();
        #endregion

        /// <summary>
        /// Given a segment, return the expected copy number - normally this is 2, but based on the reference ploidy bed file, it could be something else.  
        /// For XX samples, reference ploidy is 0 on chrY; for XY samples, reference ploidy is 1 on chrX+chrY
        /// </summary>
        public int GetReferenceCopyNumber(CanvasSegment segment)
        {
            if (!PloidyByChromosome.ContainsKey(segment.Chr)) return 2;
            int[] baseCounts = new int[5];
            baseCounts[2] = segment.End - segment.Begin;

            foreach (PloidyInterval interval in this.PloidyByChromosome[segment.Chr])
            {
                if (interval.Ploidy == 2) continue;
                int overlapStart = Math.Max(segment.Begin, interval.Start);
                if (overlapStart > segment.End) continue;
                int overlapEnd = Math.Min(segment.End, interval.End);
                int overlapBases = overlapEnd - overlapStart;
                if (overlapBases < 0) continue;
                baseCounts[2] -= overlapBases;
                baseCounts[interval.Ploidy] += overlapBases; // ASSUMPTION: Bed file ploidy shouldn't be >4 (i.e. we wouldn't handle an XXXXXY genome):
            }
            int bestCount = 0;
            int referenceCN = 2;
            for (int CN = 0; CN < baseCounts.Length; CN++)
            {
                if (baseCounts[CN] > bestCount)
                {
                    bestCount = baseCounts[CN];
                    referenceCN = CN;
                }
            }
            return referenceCN;
        }

        public static PloidyInfo LoadPloidyFromBedFile(string filePath)
        {
            PloidyInfo ploidy = new PloidyInfo();
            int count = 0;
            using (GzipReader reader = new GzipReader(filePath))
            {
                while (true)
                {
                    string fileLine = reader.ReadLine();
                    if (fileLine == null) break;
                    if (fileLine.StartsWith("##ReferenceSexChromosomeKaryotype"))
                    {
                        ploidy.HeaderLine = fileLine.Trim();
                        continue;
                    }
                    if (fileLine.Length == 0 || fileLine[0] == '#') continue;
                    string[] bits = fileLine.Split('\t');
                    string chromosome = bits[0];
                    if (!ploidy.PloidyByChromosome.ContainsKey(chromosome))
                    {
                        ploidy.PloidyByChromosome[chromosome] = new List<PloidyInterval>();
                    }
                    PloidyInterval interval = new PloidyInterval();
                    interval.Start = int.Parse(bits[1]);
                    interval.End = int.Parse(bits[2]);
                    interval.Ploidy = int.Parse(bits[4]);
                    ploidy.PloidyByChromosome[chromosome].Add(interval);
                    count++;
                }
            }
            Console.WriteLine("Reference ploidy: Loaded {0} intervals across {1} chromosomes", count, ploidy.PloidyByChromosome.Keys.Count);
            return ploidy;
        }
    }

    public class PloidyInterval
    {
        public int Start;
        public int End;
        public int Ploidy;
    }

    /// <summary>
    /// a SegmentPloidy represents a potential copy number state of a genomic interval.  It has a copy number
    /// and a major chromosome count.  
    /// </summary>
    public class SegmentPloidy
    {
        public int ID; // A 0-based index for this ploidy, for array indexing.
        public int CopyNumber;
        public int MajorChromosomeCount;
        public double MinorAlleleFrequency; // for PURE tumor
        public double MixedMinorAlleleFrequency; // for our ESTIMATED overall purity
        public double MixedCoverage; // for our ESTIMATED overall purity

        // Used by the EM algorithm
        // Dimensions: 0: MAF, 1: Coverage
        public double Omega; // weight of the Gaussian
        public double[] Mu; // means
        public double[][] Sigma; // covariance matrix

    }

    /// <summary>
    /// Represents a model point in (MAF, Coverage) space
    /// </summary>
    public class ModelPoint
    {
        public double MAF;
        public double Coverage;
        public double Weight;
        public double MAFWeight;
        public int CN;
        public int? Cluster;
        public double? KnearestNeighbour;
        public SegmentPloidy Ploidy;
        public double EmpiricalMAF;
        public double EmpiricalCoverage;
        public double Distance;
    }

    /// <summary>
    /// Represents a point in (MAF, Coverage) space to be clustered:
    /// </summary>
    public class SegmentInfo : ModelPoint
    {
        public CanvasSegment Segment;
        public Dictionary<ModelPoint, double> PosteriorProbs;

    }

    public class CoverageModel
    {
        public double Deviation;
        public double PrecisionDeviation;
        public double AccuracyDeviation;
        public double DiploidCoverage;
        public double Ploidy;
        public List<int> CNs = new List<int>();
    }

    /// <summary>
    /// Stores information about clustering of CanvasSegment objects
    /// </summary>
    public class ClusterModel
    {
        public int? ClusterID;
        public bool IsHeterogeneous = false;
        public double ClusterMedianDistance;
    }
}
