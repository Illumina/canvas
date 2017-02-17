using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Isas.SequencingFiles;
using MathNet.Numerics.Statistics;

namespace CanvasCommon
{
    public class CanvasSegmentWriter
    {
        /// <summary>
        /// Integrity check, to ensure that our reference FASTA file is in sync with our inputs.  
        /// </summary>
        private static void SanityCheckChromosomeNames(GenomeMetadata genome, List<CanvasSegment> segments)
        {
            var chromosomeNames = new HashSet<string>();
            foreach (GenomeMetadata.SequenceMetadata chromosome in genome.Sequences)
            {
                chromosomeNames.Add(chromosome.Name.ToLowerInvariant());
            }
            foreach (
                CanvasSegment segment in
                segments.Where(segment => !chromosomeNames.Contains(segment.Chr.ToLowerInvariant())))
            {
                throw new Exception($"Integrity check error: Segment found at unknown chromosome '{segment.Chr}'");
            }
        }

        public static void AddPloidyAndCoverageHeaders(BgzipOrStreamWriter writer, List<CanvasSegment> segments, double? diploidCoverage)
        {
            double totalPloidy = 0;
            double totalWeight = 0;
            foreach (CanvasSegment segment in segments.Where(segment => segment.Filter == "PASS"))
            {
                totalWeight += segment.End - segment.Begin;
                totalPloidy += segment.CopyNumber * (segment.End - segment.Begin);
            }
            if (totalWeight > 0)
            {
                writer.WriteLine($"##OverallPloidy={totalPloidy / totalWeight:F2}");
                if (diploidCoverage != null) writer.WriteLine($"##DiploidCoverage={diploidCoverage:F2}");
            }
        }

        private static GenomeMetadata WriteVcfHeader(List<CanvasSegment> segments, double? diploidCoverage,
            string wholeGenomeFastaDirectory, List<string> sampleNames, List<string> extraHeaders, int qualityThreshold,
            BgzipOrStreamWriter writer, out string denovoQualityFilter, int? denovoQualityThreshold = null)
        {
            // Write the VCF header:
            writer.WriteLine("##fileformat=VCFv4.1");
            writer.WriteLine($"##source={CanvasVersionInfo.NameString} {CanvasVersionInfo.VersionString}");
            writer.WriteLine($"##reference={Path.Combine(wholeGenomeFastaDirectory, "genome.fa")}");
            AddPloidyAndCoverageHeaders(writer, segments, diploidCoverage);
            foreach (string header in extraHeaders ?? new List<string>())
            {
                writer.WriteLine(header);
            }

            GenomeMetadata genome = new GenomeMetadata();
            genome.Deserialize(Path.Combine(wholeGenomeFastaDirectory, "GenomeSize.xml"));
            foreach (GenomeMetadata.SequenceMetadata chromosome in genome.Sequences)
            {
                writer.WriteLine($"##contig=<ID={chromosome.Name},length={chromosome.Length}>");
            }
            string qualityFilter = $"q{qualityThreshold}";
            writer.WriteLine("##ALT=<ID=CNV,Description=\"Copy number variable region\">");
            writer.WriteLine($"##FILTER=<ID={qualityFilter},Description=\"Quality below {qualityThreshold}\">");
            writer.WriteLine("##FILTER=<ID=L10kb,Description=\"Length shorter than 10kb\">");
            denovoQualityFilter = "";
            if (denovoQualityThreshold.HasValue)
            {
                denovoQualityFilter = $"dq{denovoQualityThreshold}";
                writer.WriteLine(
                    $"##FILTER=<ID={denovoQualityFilter},Description=\"De novo quality score above {denovoQualityThreshold.Value}\">");
            }
            writer.WriteLine(
                "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">");
            writer.WriteLine(
                "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">");
            writer.WriteLine(
                "##INFO=<ID=CNVLEN,Number=1,Type=Integer,Description=\"Number of reference positions spanned by this CNV\">");
            writer.WriteLine(
                "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">");
            writer.WriteLine("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">");
            writer.WriteLine(
                "##INFO=<ID=DQSCORE,Number=1,Type=String,Description=\"De novo Phred-scaled quality score\">");
            writer.WriteLine("##INFO=<ID=SUBCLONAL,Number=0,Type=Flag,Description=\"Subclonal variant\">");
            writer.WriteLine("##FORMAT=<ID=RC,Number=1,Type=Float,Description=\"Mean counts per bin in the region\">");
            writer.WriteLine("##FORMAT=<ID=BC,Number=1,Type=Float,Description=\"Number of bins in the region\">");
            writer.WriteLine(
                "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number genotype for imprecise events\">");
            writer.WriteLine(
                "##FORMAT=<ID=MCC,Number=1,Type=Integer,Description=\"Major chromosome count (equal to copy number for LOH regions)\">");
            string names = sampleNames.Count == 1 ? sampleNames.First() : string.Join("\t", sampleNames.ToArray()) ;
            writer.WriteLine("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + names);
            SanityCheckChromosomeNames(genome, segments);
            return genome;
        }


        /// <summary>
        /// Outputs the copy number calls to a text file.
        /// </summary>
        private static void WriteVariants(IReadOnlyCollection<List<CanvasSegment>> segments, PloidyInfo ploidy, GenomeMetadata genome,
            BgzipOrStreamWriter writer, string denovoQualityFilter, int? denovoQualityThreshold = null)
        {
            var nSamples = segments.Count;
            foreach (GenomeMetadata.SequenceMetadata chromosome in genome.Sequences)
            {
                for (int segmentIndex = 0; segmentIndex < segments.Single().Count; segmentIndex++)
                {
                    if (!segments.First()[segmentIndex].Chr.Equals(chromosome.Name, StringComparison.OrdinalIgnoreCase))
                        continue;
                    var referenceCopyNumbers = segments.Select(segment => ploidy?.GetReferenceCopyNumber(segment[segmentIndex]) ?? 2).ToList();
                    var currentSegments = segments.Select(x => x[segmentIndex]).ToList();
                    var cnvTypes = new List<CnvType>();
                    for (int sampleIndex = 0; sampleIndex < nSamples; sampleIndex++)
                    {
                        cnvTypes.Add(currentSegments[sampleIndex].GetCnvType(referenceCopyNumbers[sampleIndex]));
                    }
                    CnvType cnvType;
                    if (cnvTypes.TrueForAll(x => x == CnvType.Reference))
                        cnvType = CnvType.Reference;
                    else if (cnvTypes.TrueForAll(x => x == CnvType.Reference | x == CnvType.Loss))
                        cnvType = CnvType.Loss;
                    else if (cnvTypes.TrueForAll(x => x == CnvType.Reference | x == CnvType.Gain))
                        cnvType = CnvType.Gain;
                    else if (cnvTypes.TrueForAll(x => x == CnvType.Reference | x == CnvType.LossOfHeterozygosity))
                        cnvType = CnvType.LossOfHeterozygosity;
                    else
                        cnvType = CnvType.ComplexCnv;
                            
                    WriteVcfVariantInfo(ploidy, writer, segments.First()[segmentIndex], cnvType);
                    //  FORMAT field
                    if (segments.Count == 1)
                        WriteSingleSampleInfo(writer, segments.First()[segmentIndex]);
                    else
                        WriteMultiSampleInfo(writer, segments.Select(x=>x[segmentIndex]).ToList());
                }
            }
        }

        private static void WriteSingleSampleInfo(BgzipOrStreamWriter writer, CanvasSegment segment)
        {
            writer.Write("\tRC:BC:CN", segment.End);
            if (segment.MajorChromosomeCount.HasValue)
            {
                writer.Write(":MCC");
            }
            writer.Write("\t{1}:{2}:{3}", segment.End, Math.Round(segment.MeanCount, 0, MidpointRounding.AwayFromZero), 
                segment.BinCount, segment.CopyNumber);
            if (segment.MajorChromosomeCount.HasValue)
            {
                writer.Write(":{0}", segment.MajorChromosomeCount);
            }
            writer.WriteLine();
        }

        private static void WriteMultiSampleInfo(BgzipOrStreamWriter writer, List<CanvasSegment> segments)
        {
            writer.Write("\tRC:BC:CN:MCC:DQSCORE");
            string nullValue = ".";
            foreach (var segment in segments)
            {
                var mcc = segment.MajorChromosomeCount.HasValue ? segment.MajorChromosomeCount.ToString() : nullValue;
                var dqscore = segment.DQScore.HasValue ? $"{segment.DQScore.Value:F2}" : nullValue;
                var rc = Math.Round(segment.MeanCount, 0, MidpointRounding.AwayFromZero);
                writer.Write($"{rc}:{segment.BinCount}:{ segment.CopyNumber}:{mcc}:{dqscore}");
            }
        }

        /// <summary>
        /// Write to a file a single CanvasSegment record as a non-Format VCF columns 
        /// </summary>
        /// <param name="ploidy"></param>
        /// <param name="writer"></param>
        /// <param name="denovoQualityFilter"></param>
        /// <param name="denovoQualityThreshold"></param>
        /// <param name="segment"></param>
        /// <param name="chromosome"></param>
        /// <returns></returns>
        private static void WriteVcfVariantInfo(PloidyInfo ploidy, BgzipOrStreamWriter writer, CanvasSegment segment, CnvType cnvType)
        {

            // From vcf 4.1 spec:
            //     If any of the ALT alleles is a symbolic allele (an angle-bracketed ID String “<ID>”) then the padding base is required and POS denotes the 
            //     coordinate of the base preceding the polymorphism.
            string alternateAllele = cnvType.ToAltId();
            int position = (alternateAllele.StartsWith("<") && alternateAllele.EndsWith(">"))
                ? segment.Begin
                : segment.Begin + 1;

            writer.Write($"{segment.Chr}\t{position}\tCanvas:{cnvType.ToVcfId()}:{segment.Chr}:{segment.Begin + 1}-{segment.End}\t");
            writer.Write($"N\t{alternateAllele}\t{segment.QScore:F2}\t{segment.Filter}\t", alternateAllele, segment.QScore, segment.Filter);

            if (cnvType != CnvType.Reference)
                writer.Write($"SVTYPE={cnvType.ToSvType()};");
            if (segment.IsHeterogeneous)
                writer.Write("SUBCLONAL;");

            writer.Write($"END={segment.End}");
            if (cnvType != CnvType.Reference)
                writer.Write($";CNVLEN={segment.End - segment.Begin}");

            if (segment.StartConfidenceInterval != null)
                writer.Write($";CIPOS={segment.StartConfidenceInterval.Item1},{segment.StartConfidenceInterval.Item2}");

            if (segment.EndConfidenceInterval != null)
                writer.Write($";CIEND={segment.EndConfidenceInterval.Item1},{segment.EndConfidenceInterval.Item2}");
        }


        public static void WriteSegments(string outVcfPath, List<CanvasSegment> segments, double? diploidCoverage,
                string wholeGenomeFastaDirectory, string sampleName,
                List<string> extraHeaders, PloidyInfo ploidy, int qualityThreshold, int? denovoQualityThreshold = null)
        {
            using (BgzipOrStreamWriter writer = new BgzipOrStreamWriter(outVcfPath))
            {
                string denovoQualityFilter;
                var genome = WriteVcfHeader(segments, diploidCoverage, wholeGenomeFastaDirectory, new List<string> {sampleName},
                    extraHeaders, qualityThreshold, writer, out denovoQualityFilter, denovoQualityThreshold);
                WriteVariants(new List<List<CanvasSegment>> {segments.ToList()}, ploidy, genome, writer, denovoQualityFilter, denovoQualityThreshold);
            }
        }
        public static void WriteMultiSampleSegments(string outVcfPath, List<List<CanvasSegment>> segments, List<double?> diploidCoverage,
        string wholeGenomeFastaDirectory, List<string> sampleNames,
        List<string> extraHeaders, PloidyInfo ploidy, int qualityThreshold, int? denovoQualityThreshold = null)
        {
            using (BgzipOrStreamWriter writer = new BgzipOrStreamWriter(outVcfPath))
            {
                string denovoQualityFilter;
                var genome = WriteVcfHeader(segments.First(), diploidCoverage.Mean(), wholeGenomeFastaDirectory, sampleNames,
                    extraHeaders, qualityThreshold, writer, out denovoQualityFilter, denovoQualityThreshold);
                WriteVariants(segments, ploidy, genome, writer, denovoQualityFilter, denovoQualityThreshold);
            }
        }
    }
}