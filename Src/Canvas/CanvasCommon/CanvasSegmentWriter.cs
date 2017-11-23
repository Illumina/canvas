using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Isas.Framework.DataTypes;
using Isas.Framework.DataTypes.Maps;
using Isas.SequencingFiles;


namespace CanvasCommon
{
    public class CanvasSegmentWriter
    {
        /// <summary>
        /// Integrity check, to ensure that our reference FASTA file is in sync with our inputs.  
        /// </summary>
        private static void SanityCheckChromosomeNames(GenomeMetadata genome, IEnumerable<CanvasSegment> segments)
        {
            var chromosomeNames = new HashSet<string>();
            foreach (GenomeMetadata.SequenceMetadata chromosome in genome.Contigs())
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
            foreach (CanvasSegment segment in segments.Where(segment => segment.Filter.IsPass))
            {
                totalWeight += segment.Length;
                totalPloidy += segment.CopyNumber * (segment.Length);
            }
            if (totalWeight > 0)
            {
                writer.WriteLine($"##OverallPloidy={totalPloidy / totalWeight:F2}");
                if (diploidCoverage != null) writer.WriteLine($"##DiploidCoverage={diploidCoverage:F2}");
            }
        }

        private static GenomeMetadata WriteVcfHeader(List<CanvasSegment> segments, double? diploidCoverage,
            string wholeGenomeFastaDirectory, List<string> sampleNames, List<string> extraHeaders, BgzipOrStreamWriter writer, int qualityThreshold,
            int? denovoQualityThreshold, int? sizeThreshold)
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
            genome.Deserialize(new FileLocation(Path.Combine(wholeGenomeFastaDirectory, "GenomeSize.xml")));

            foreach (GenomeMetadata.SequenceMetadata chromosome in genome.Contigs())
            {
                writer.WriteLine($"##contig=<ID={chromosome.Name},length={chromosome.Length}>");
            }
            string qualityFilter = $"q{qualityThreshold}";
            writer.WriteLine("##ALT=<ID=CNV,Description=\"Copy number variable region\">");
            WriteHeaderAllAltCnTags(writer);
            writer.WriteLine($"##FILTER=<ID={qualityFilter},Description=\"Quality below {qualityThreshold}\">");
            if (sizeThreshold.HasValue)
            {
                string sizeFilterValue = CanvasFilter.FormatCnvSizeWithSuffix(sizeThreshold.Value);
                writer.WriteLine($"##FILTER=<ID=L{sizeFilterValue},Description=\"Length shorter than {sizeFilterValue}\">");
            }
            writer.WriteLine("##FILTER=<ID=FailedFT,Description=\"Sample-level filter failed in all the samples\">");
            writer.WriteLine("##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">");
            writer.WriteLine("##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">");
            writer.WriteLine("##INFO=<ID=CNVLEN,Number=1,Type=Integer,Description=\"Number of reference positions spanned by this CNV\">");
            writer.WriteLine("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">");
            writer.WriteLine("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">");
            writer.WriteLine("##INFO=<ID=SUBCLONAL,Number=0,Type=Flag,Description=\"Subclonal variant\">");
            writer.WriteLine("##INFO=<ID=COMMONCNV,Number=0,Type=Flag,Description=\"Common CNV variant identified from pre-specified bed intervals\">");
            writer.WriteLine("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
            writer.WriteLine("##FORMAT=<ID=RC,Number=1,Type=Float,Description=\"Mean counts per bin in the region\">");
            writer.WriteLine("##FORMAT=<ID=BC,Number=1,Type=Float,Description=\"Number of bins in the region\">");
            writer.WriteLine("##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number genotype for imprecise events\">");
            writer.WriteLine("##FORMAT=<ID=MCC,Number=1,Type=Integer,Description=\"Major chromosome count (equal to copy number for LOH regions)\">");
            writer.WriteLine("##FORMAT=<ID=MCCQ,Number=1,Type=Float,Description=\"Major chromosome count quality score\">");
            writer.WriteLine("##FORMAT=<ID=QS,Number=1,Type=Float,Description=\"Phred-scaled quality score. If CN is reference then this is -10log10(prob(variant)) otherwise this is -10log10(prob(no variant).\">");
            if (denovoQualityThreshold.HasValue)
            {
                writer.WriteLine($"##FORMAT=<ID=DQ,Number=1,Type=Float,Description=\"De novo quality. Threshold for passing de novo call: {denovoQualityThreshold}\">");
            }
            writer.WriteLine("##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Sample-level filter\">");
            var titleColumns = new List<string> { "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT" };
            titleColumns.AddRange(sampleNames);
            writer.WriteLine(string.Join("\t", titleColumns));
            SanityCheckChromosomeNames(genome, segments);
            return genome;
        }

        public static void WriteHeaderAllAltCnTags(BgzipOrStreamWriter writer, int maxCopyNum = 5)
        {
            foreach (var copyNum in Enumerable.Range(0, maxCopyNum + 1))
            {
                if (copyNum == 1) continue;
                writer.WriteLine($"##ALT=<ID=CN{copyNum},Description=\"Copy number allele: {copyNum} copies\">");
            }
        }

        /// <summary>
        /// Outputs the copy number calls to a text file.
        /// </summary>
        private static void WriteVariants(IEnumerable<ISampleMap<CanvasSegment>> segmentsOfAllSamples, List<PloidyInfo> ploidies, GenomeMetadata genome,
            BgzipOrStreamWriter writer, bool isPedigreeInfoSupplied = true, int? denovoQualityThreshold = null)
        {
            var segmentsOfAllSamplesArray = segmentsOfAllSamples.ToArray(); // TODO: not necessary when chrom match logic has been updated
            int nSamples = segmentsOfAllSamplesArray.First().Values.Count();
            foreach (GenomeMetadata.SequenceMetadata chromosome in genome.Contigs()) //TODO: this is extremely inefficient. Segments should be sorted by chromosome
            {
                foreach (var sampleMap in segmentsOfAllSamplesArray)
                {
                    var currentSegments = sampleMap.Values.ToArray();
                    var firstSampleSegment = currentSegments.First();
                    if (!firstSampleSegment.Chr.Equals(chromosome.Name, StringComparison.OrdinalIgnoreCase)
                    ) //TODO: this is extremely inefficient. Segments should be sorted by chromosome
                        continue;
                    var recordLevelFilter = CanvasFilter.GetRecordLevelFilterFromSampleFiltersOnly(
                                            sampleMap
                                            .Select(x => x.Value.Filter)
                                            .ToReadOnlyList())
                                            .ToVcfString();
                    var referenceCopyNumbers = currentSegments.Zip(ploidies,
                        (segment, ploidy) => ploidy?.GetReferenceCopyNumber(segment) ?? 2).ToList();
                    var cnvTypes = new CnvType[nSamples];
                    var sampleSetAlleleCopyNumbers = new int[nSamples][];
                    for (int sampleIndex = 0; sampleIndex < nSamples; sampleIndex++)
                    {
                        (cnvTypes[sampleIndex], sampleSetAlleleCopyNumbers[sampleIndex]) = currentSegments[sampleIndex]
                            .GetCnvTypeAndAlleleCopyNumbers(referenceCopyNumbers[sampleIndex]);
                    }
                    var sampleSetCnvType = AssignCnvType(cnvTypes);
                    var (alternateAllele, genotypes) = GetAltAllelesAndGenotypes(sampleSetAlleleCopyNumbers);
                    WriteColumnsUntilInfoField(writer, firstSampleSegment, sampleSetCnvType, alternateAllele,
                        recordLevelFilter, nSamples > 1);
                    WriteFormatAndSampleFields(writer, currentSegments, genotypes,
                        denovoQualityThreshold.HasValue);
                }
            }
        }

        private static CnvType AssignCnvType(CnvType[] cnvTypes)
        {
            var nonRefTypes = cnvTypes.Where(x => x != CnvType.Reference).Distinct().ToArray();
            if (!nonRefTypes.Any()) return CnvType.Reference;
            if (nonRefTypes.Length > 1) return CnvType.ComplexCnv;
            return nonRefTypes.First();
        }

        internal static (string, string[]) GetAltAllelesAndGenotypes(int[][] sampleSetAlleleCopyNumbers)
        {
            var uniqAltAlleles = sampleSetAlleleCopyNumbers.SelectMany(x => x).Distinct().Where(x => x != 1 && x != -1).OrderBy(x => x).ToList();
            var altAlleles = uniqAltAlleles.Select(x => $"<CN{x}>").ToArray();
            string altAlleleString = ".";
            if (altAlleles.Any())
            {
                if (altAlleles.Last() == $"<CN{int.MaxValue}>")
                    altAlleles[altAlleles.Length - 1] = "<DUP>"; // DUP is always the last one
                altAlleleString = string.Join(",", altAlleles);
            }
            var sampleGenotypes = sampleSetAlleleCopyNumbers.Select(x => AlleleCopyNumberToGenotype(x, uniqAltAlleles)).ToArray();
            return (altAlleleString, sampleGenotypes);
        }


        private static string AlleleCopyNumberToGenotype(int[] alleleCopyNumbers, List<int> uniqAltAlleles)
        {
            var genotypes = new string[alleleCopyNumbers.Length];
            for (int i = 0; i < alleleCopyNumbers.Length; i++)
            {
                switch (alleleCopyNumbers[i])
                {
                    case 1:
                        genotypes[i] = "0";
                        break;
                    case -1:
                        genotypes[i] = ".";
                        break;
                    default:
                        genotypes[i] = (uniqAltAlleles.IndexOf(alleleCopyNumbers[i]) + 1).ToString();
                        break;
                }
            }
            return string.Join("/", genotypes.OrderBy(GenotypeToIntMapping));
        }

        private static int GenotypeToIntMapping(string genotype)  // '.' first, then in numeric order
        {
            return genotype == "." ? -1 : int.Parse(genotype);
        }

        private static void WriteFormatAndSampleFields(BgzipOrStreamWriter writer, CanvasSegment[] segments, string[] genotypes, bool reportDQ)
        {
            const string nullValue = ".";
            writer.Write("\tGT:RC:BC:CN:MCC:MCCQ:QS:FT");
            if (reportDQ)
                writer.Write(":DQ");
            for (int i = 0; i < segments.Length; i++)
            {
                var segment = segments[i];
                string mcc = segment.MajorChromosomeCount.HasValue ? segment.MajorChromosomeCount.ToString() : nullValue;
                string mccq = segment.MajorChromosomeCountScore.HasValue ? $"{segment.MajorChromosomeCountScore.Value:F2}" : nullValue;
       
                writer.Write($"\t{genotypes[i]}:{segment.MedianCount:F2}:{segment.BinCount}:{segment.CopyNumber}:{mcc}:{mccq}:{segment.QScore:F2}:{segment.Filter.ToVcfString()}");
                if (reportDQ)
                {
                    string dqscore = segment.DqScore.HasValue ? $"{segment.DqScore.Value:F2}" : nullValue;
                    writer.Write($":{dqscore}");
                }
                writer.WriteLine();
            }
        }


        /// <summary>
        /// Write to a file a single CanvasSegment record as a non-sample VCF columns 
        /// </summary>
        /// <param name="writer"></param>
        /// <param name="firstSampleSegment"></param>
        /// <param name="alternateAllele"></param>
        /// <param name="recordLevelFilter"></param>
        /// <param name="sampleSetCnvType"></param>
        /// <param name="isMultisample"></param>
        /// <returns></returns>
        private static void WriteColumnsUntilInfoField(BgzipOrStreamWriter writer, CanvasSegment firstSampleSegment, CnvType sampleSetCnvType, string alternateAllele, string recordLevelFilter, bool isMultisample)
        {
            // From vcf 4.1 spec:
            //     If any of the ALT alleles is a symbolic allele (an angle-bracketed ID String “<ID>”) then the padding base is required and POS denotes the 
            //     coordinate of the base preceding the polymorphism.
            int position = (alternateAllele.StartsWith("<") && alternateAllele.EndsWith(">"))
                ? firstSampleSegment.Begin
                : firstSampleSegment.Begin + 1;
            writer.Write($"{firstSampleSegment.Chr}\t{position}\tCanvas:{sampleSetCnvType.ToVcfId()}:{firstSampleSegment.Chr}:{firstSampleSegment.Begin + 1}-{firstSampleSegment.End}\t");
            string qScore = isMultisample ? "." : $"{firstSampleSegment.QScore:F2}";
            writer.Write($"N\t{alternateAllele}\t{qScore}\t{recordLevelFilter}\t");

            if (sampleSetCnvType != CnvType.Reference)
                writer.Write($"SVTYPE={sampleSetCnvType.ToSvType()};");

            if (firstSampleSegment.IsHeterogeneous)
                writer.Write("SUBCLONAL;");

            if (firstSampleSegment.IsCommonCnv)
                writer.Write("COMMONCNV;");

            writer.Write($"END={firstSampleSegment.End}");

            if (sampleSetCnvType != CnvType.Reference)
                writer.Write($";CNVLEN={firstSampleSegment.Length}");

            if (firstSampleSegment.StartConfidenceInterval != null)
                writer.Write($";CIPOS={firstSampleSegment.StartConfidenceInterval.Item1},{firstSampleSegment.StartConfidenceInterval.Item2}");
            if (firstSampleSegment.EndConfidenceInterval != null)
                writer.Write($";CIEND={firstSampleSegment.EndConfidenceInterval.Item1},{firstSampleSegment.EndConfidenceInterval.Item2}");
        }


        public static void WriteSegments(string outVcfPath, List<CanvasSegment> segments, double? diploidCoverage,
                string wholeGenomeFastaDirectory, string sampleName,
                List<string> extraHeaders, PloidyInfo ploidy, int qualityThreshold, bool isPedigreeInfoSupplied, int? denovoQualityThreshold, int? sizeThreshold)
        {
            using (BgzipOrStreamWriter writer = new BgzipOrStreamWriter(outVcfPath))
            {
                var genome = WriteVcfHeader(segments, diploidCoverage, wholeGenomeFastaDirectory, new List<string> { sampleName },
                    extraHeaders, writer, qualityThreshold, denovoQualityThreshold, sizeThreshold);
                var sampleId = new SampleId(sampleName);
                var segmentsOfAllSamples = segments.Select(x => new SampleMap<CanvasSegment> {{sampleId, x}});
                WriteVariants(segmentsOfAllSamples, new List<PloidyInfo> { ploidy }, genome, writer, isPedigreeInfoSupplied, denovoQualityThreshold);
            }
        }

        public static void WriteMultiSampleSegments(string outVcfPath, ISampleMap<List<CanvasSegment>> segments, List<double> diploidCoverage,
        string wholeGenomeFastaDirectory, List<string> sampleNames, List<string> extraHeaders, List<PloidyInfo> ploidies, int qualityThreshold, int? denovoQualityThreshold, int? sizeThreshold, bool isPedigreeInfoSupplied = true)
        {
            using (BgzipOrStreamWriter writer = new BgzipOrStreamWriter(outVcfPath))
            {
                var genome = WriteVcfHeader(segments.Values.First(), diploidCoverage.Average(), wholeGenomeFastaDirectory, sampleNames,
                    extraHeaders, writer, qualityThreshold, denovoQualityThreshold, sizeThreshold);
                WriteVariants(segments.Zip(), ploidies, genome, writer, isPedigreeInfoSupplied, denovoQualityThreshold);
            }
        }
    }
}
