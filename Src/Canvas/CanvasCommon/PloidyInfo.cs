using System;
using System.Collections.Generic;
using System.Linq;
using Illumina.Common;
using Isas.SequencingFiles;
using Isas.SequencingFiles.Vcf;

namespace CanvasCommon
{
    public class PloidyInfo
    {
        #region Members

        public const int OutlierClusterFlag = -1;
        public const int UndersegmentedClusterFlag = -2;
        public string HeaderLine;
        public Dictionary<string, List<PloidyInterval>> PloidyByChromosome = new Dictionary<string, List<PloidyInterval>>();

        /// <summary>
        /// Make sure the list of ploidy intervals for a particular chromosome in the PloidyByChromosome dictionary 
        /// can be found for either chromosome naming convention (e.g. "chrX" and "X") 
        /// Also make sure that each input chromosome appears in the PloidyByChromosome dictionary.
        /// If there is currently no entry for a particular chromosome add an empty List of ploidy intervals 
        /// </summary>
        public void MakeChromsomeNameAgnosticWithAllChromosomes(IEnumerable<string> chromosomes)
        {
            Dictionary<string, List<PloidyInterval>> ploidyByChromosomeAllChromosomes = new Dictionary<string, List<PloidyInterval>>(PloidyByChromosome);
            foreach (var chromosome in chromosomes)
            {
                var altChromosome = chromosome.StartsWith("chr") ? chromosome.Substring(3) : "chr" + chromosome;
                var ploidyIntervals = new List<PloidyInterval>();
                if (ploidyByChromosomeAllChromosomes.ContainsKey(chromosome))
                {
                    ploidyIntervals = ploidyByChromosomeAllChromosomes[chromosome];
                }
                else if (ploidyByChromosomeAllChromosomes.ContainsKey(altChromosome))
                {
                    ploidyIntervals = ploidyByChromosomeAllChromosomes[altChromosome];
                }
                if (!ploidyByChromosomeAllChromosomes.ContainsKey(chromosome))
                {
                    ploidyByChromosomeAllChromosomes.Add(chromosome, ploidyIntervals);
                }
                if (!ploidyByChromosomeAllChromosomes.ContainsKey(altChromosome))
                {
                    ploidyByChromosomeAllChromosomes.Add(altChromosome, ploidyIntervals);
                }
            }
            PloidyByChromosome = ploidyByChromosomeAllChromosomes;
        }

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
                baseCounts[interval.Ploidy] += overlapBases; // ASSUMPTION: Vcf file ploidy shouldn't be >4 (i.e. we wouldn't handle an XXXXXY genome):
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

        // Only one sample column allowed, if no sampleId provided, 
        public static PloidyInfo LoadPloidyFromVcfFileNoSampleId(string vcfPath)
        {
            // check how many sample columns in the VCF file
            using (VcfReader reader = new VcfReader(vcfPath))
            {
                var sampleCount = reader.Samples.Count;
                if (sampleCount == 0) 
                    throw new ArgumentException(
                        $"File '{vcfPath}' does not contain any genotype column");
                else if (sampleCount > 1)
                    throw new ArgumentException(
                        $"File '{vcfPath}' cannot have more than one genotype columns when no sample ID provided'");
            }
            return LoadPloidyFromVcfFile(vcfPath, 0);
        }

        private static PloidyInfo LoadPloidyFromVcfFile(string vcfPath, int sampleIndex) 
        {
            PloidyInfo ploidy = new PloidyInfo();

            using (VcfReader reader = new VcfReader(vcfPath))
            {
                ploidy.HeaderLine = string.Join(" ", reader.HeaderLines);

                while (true)
                {
                    VcfVariant record;
                    bool result = reader.GetNextVariant(out record);
                    if (!result) break;
                    string chromosome = record.ReferenceName;
                    if (!ploidy.PloidyByChromosome.ContainsKey(chromosome))
                    {
                        ploidy.PloidyByChromosome[chromosome] = new List<PloidyInterval>();
                    }
                    PloidyInterval interval = new PloidyInterval(chromosome);
                    interval.Start = record.ReferencePosition;
                    interval.End = int.Parse(record.InfoFields["END"]);
                    var genotypeColumn = record.GenotypeColumns[sampleIndex];
                    if (genotypeColumn.ContainsKey("CN"))
                    {
                        var value = genotypeColumn["CN"];
                        interval.Ploidy = value == "." ? 2 : int.Parse(value);
                    }
                    else
                        throw new ArgumentException($"File '{vcfPath}' must contain one genotype CN column!");
                    ploidy.PloidyByChromosome[chromosome].Add(interval);
                }
            }
            return ploidy;
        }


        public static PloidyInfo LoadPloidyFromVcfFile(string vcfPath, string sampleId)
        {
            int sampleIndex;
            using (VcfReader reader = new VcfReader(vcfPath))
            {
                sampleIndex = reader.Samples.IndexOf(sampleId);
                if (sampleIndex == -1)
                    throw new ArgumentException(
                        $"File '{vcfPath}' does not contain a genotype column for sample '{sampleId}'");
            }
            return LoadPloidyFromVcfFile(vcfPath, sampleIndex);
        }
    }

    public class PloidyInterval
    {
        public string Chromosome { get; }
        public int Start;
        public int End;
        public int Ploidy;

        public PloidyInterval(string chromosome)
        {
            Chromosome = chromosome;
        }

        public override string ToString()
        {
            return $"{Chromosome}:{Start + 1}-{End}";
        }
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
        public int? ClusterId;
        public int? FinalClusterId;
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
        public ClusterInfo Cluster;
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

}
