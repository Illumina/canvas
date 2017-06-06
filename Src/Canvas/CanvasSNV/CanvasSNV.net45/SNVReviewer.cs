using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Illumina.Common;
using Illumina.Common.CSV;
using Isas.SequencingFiles;
using Isas.SequencingFiles.Vcf;

namespace CanvasSNV
{
    /// <summary>
    /// The SNV reviewer takes as input a .vcf file (from the normal sample) and a .bam file (from the paired tumor sample).
    /// It processes a particular chromosome at a time.
    /// For all variants in the .vcf file that pass our filters, we review the variant allele frequency in the .bam file
    /// from the tumor sample.  We output that data in tabular format.
    /// </summary>
    public class SNVReviewer
    {
        #region Members
        protected readonly string Chromosome;
        protected readonly string VcfPath;
        protected readonly string BamPath;
        protected readonly string OutputPath;
        protected readonly string SampleName;
        protected readonly bool IsDbSnpVcf;
        private readonly bool IsSomatic;
        List<VcfVariant> Variants;
        int[] ReferenceCounts;
        int[] VariantCounts;
        int MinimumBaseQScore = 20; // You must be >= this base call quality to be counted.
        protected readonly int MinimumMapQ;

        #endregion

        public SNVReviewer(string chromosome, string vcfPath, string bamPath, string outputPath, string sampleName, bool isDbSnpVcf, int minMapQ, bool isSomatic)
        {
            SampleName = sampleName;
            Chromosome = chromosome;
            VcfPath = vcfPath;
            BamPath = bamPath;
            OutputPath = outputPath;
            IsDbSnpVcf = isDbSnpVcf;
            MinimumMapQ = minMapQ;
            IsSomatic = isSomatic;
        }

        public int Run()
        {

            if (!File.Exists(VcfPath))
            {
                Console.Error.WriteLine("Error: Input vcf file not found at {0}", VcfPath);
                return 1;
            }
            if (!File.Exists(BamPath))
            {
                Console.Error.WriteLine("Error: Input bam file not found at {0}", BamPath);
                return 1;
            }

            this.LoadVariants(VcfPath, IsSomatic);
            ReferenceCounts = new int[Variants.Count];
            VariantCounts = new int[Variants.Count];
            this.ProcessBamFile(BamPath);
            this.WriteResults(OutputPath);
            return 0;
        }

        /// <summary>
        /// Step 1: Load the normal het SNVs of interest.
        /// </summary>
        protected void LoadVariants(string vcfPath, bool isSomatic)
        {
            Console.WriteLine("{0} Loading variants of interest from {1}", DateTime.Now, vcfPath);
            this.Variants = new List<VcfVariant>();
            int overallCount = 0;
            int countThisChromosome = 0;
            int sampleIndex = 0;
            using (VcfReader reader = new VcfReader(vcfPath, requireGenotypes: false))
            {
                if (!SampleName.IsNullOrEmpty() && IsDbSnpVcf == false)
                {
                    if (reader.Samples.All(sample => sample != SampleName))
                        throw new ArgumentException($"File '{vcfPath}' should contain one genotypes column corresponding to sample {SampleName}");
                    sampleIndex = reader.Samples.IndexOf(SampleName);
                }
                else
                {
                    if (reader.Samples.Count > 1)
                        throw new ArgumentException($"File '{vcfPath}' conatins >1 samples, name for a sample of interest must be provided");
                }

                VcfVariant variant = new VcfVariant();
                while (true)
                {
                    bool result = reader.GetNextVariant(out variant);
                    if (!result) break;
                    overallCount++;
                    if (variant.ReferenceName != this.Chromosome)
                    {
                        // Shortcut: If we've seen records for the desired chromosome, then as soon as we hit another chromosome,
                        // we can abort:
                        if (countThisChromosome > 0) break;
                        continue;
                    }
                    countThisChromosome++;
                    // Single-allele SNVs only:
                    if (variant.VariantAlleles.Length != 1 || variant.VariantAlleles[0].Length != 1 || variant.ReferenceAllele.Length != 1) continue;
                    // PF variants only:
                    if (variant.GenotypeColumns != null && variant.GenotypeColumns.Any() && variant.Filters != "PASS") continue; // FILTER may not say PASS for a dbSNP VCF file
                    if (variant.GenotypeColumns != null && variant.GenotypeColumns.Any()) // not available if we use a dbSNP VCF file
                    {
                        if (!variant.GenotypeColumns[sampleIndex].ContainsKey("GT")) continue; // no genotype - we don't know if it's a het SNV.
                        if (isSomatic)
                        {
                            string genotype = variant.GenotypeColumns[0]["GT"];
                            if (genotype != "0/1" && genotype != "1/0" && genotype != "0|1" && genotype != "1|0")
                                continue;
                        }

                        // Also require they have a high enough quality score:
                        if (variant.GenotypeColumns[sampleIndex].ContainsKey("GQX")) // Note: Allow no GQX field, in case we want to use another caller (e.g. Pisces) and not crash
                        {
                            if (variant.GenotypeColumns[sampleIndex]["GQX"].Equals("."))
                                continue;
                            if (float.Parse(variant.GenotypeColumns[sampleIndex]["GQX"]) < 30)
                                continue;
                        }
                    }
                    // Note: Let's NOT require the variant be in dbSNP.  Maybe we didn't do annotation, either because
                    // we chose not to or because we're on a reference without annotation available.
                    //if (variant.Identifier == ".") continue;
                    // Remember all the variants that pass all our tests:

                    // Empty GenotypeColumns to save space as they are no longer needed
                    variant.GenotypeColumns?.Clear();

                    this.Variants.Add(variant);
                    variant = new VcfVariant();
                }
            }
            Console.WriteLine("Retained {0} variants, out of {1} records for {2}", this.Variants.Count, countThisChromosome, this.Chromosome);
        }

        /// <summary>
        /// Step 2: Get the ref and variant allele frequencies for the variants of interest, in the tumor bam file.
        /// </summary>
        protected void ProcessBamFile(string bamPath)
        {
            Console.WriteLine("{0} Looping over bam records from {1}", DateTime.Now, bamPath);
            int overallCount = 0;
            int nextVariantIndex = 0;
            using (BamReader reader = new BamReader(bamPath))
            {
                BamAlignment read = new BamAlignment();
                int refID = reader.GetReferenceIndex(this.Chromosome);
                if (refID < 0)
                {
                    throw new ArgumentException(string.Format("Error: Chromosome name '{0}' does not match bam file at '{1}'", this.Chromosome, bamPath));
                }
                Console.WriteLine("Jump to refid {0} {1}", refID, this.Chromosome);
                reader.Jump(refID, 0);
                while (true)
                {
                    bool result = reader.GetNextAlignment(ref read, false);
                    if (!result) break;
                    if (!read.HasPosition() || read.RefID > refID) break; // We're past our chromosome of interest.
                    if (read.RefID < refID) continue; // We're not yet on our chromosome of interest.
                    overallCount++;
                    if (overallCount % 1000000 == 0)
                    {
                        Console.WriteLine("Record {0} at {1}...", overallCount, read.Position);
                    }

                    // Skip over unaligned or other non-count-worthy reads:
                    if (!read.IsPrimaryAlignment()) continue;
                    if (!read.IsMapped()) continue;
                    if (read.IsDuplicate()) continue;
                    if (read.MapQuality <= MinimumMapQ) continue;

                    // Scan forward through the variants list, to keep up with our reads:
                    while (nextVariantIndex < this.Variants.Count && this.Variants[nextVariantIndex].ReferencePosition < read.Position)
                    {
                        nextVariantIndex++;
                    }
                    if (nextVariantIndex >= this.Variants.Count) break;

                    // If the read doesn't look like it has a reasonable chance of touching the next variant, continue:
                    if (read.Position + 1000 < this.Variants[nextVariantIndex].ReferencePosition) continue;

                    // This read potentially overlaps next variant (and further variants).  Count bases!
                    ProcessReadBases(read, nextVariantIndex);
                }
            }
            Console.WriteLine("Looped over {0} bam records in all", overallCount);
        }

        /// <summary>
        /// Use the CIGAR string to map bases to chromosome positions, and check whether we see the ref base or the 
        /// variant allele for our variants of interest.
        /// </summary>
        private void ProcessReadBases(BamAlignment read, int nextVariantIndex)
        {
            int position = read.Position;
            int baseIndex = 0;
            int cigarCount = read.CigarData.Count;
            for (int opIndex = 0; opIndex < cigarCount; opIndex++)
            {
                CigarOp cigar = read.CigarData[opIndex];
                switch (cigar.Type)
                {
                    case 'M':
                        // Loop over matches/mismatches:
                        for (int index = 0; index < cigar.Length; index++, position++, baseIndex++)
                        {
                            for (int varIndex = nextVariantIndex; varIndex < this.Variants.Count; varIndex++)
                            {
                                VcfVariant variant = this.Variants[varIndex];
                                // Subtract 1: Vcf positions are 1-based, bam file positions are 0-based:
                                if (variant.ReferencePosition - 1 > position) break;
                                if (variant.ReferencePosition - 1 < position)
                                {
                                    nextVariantIndex++;
                                    continue;
                                }
                                if (read.Qualities[baseIndex] < MinimumBaseQScore) continue; // Skip low-quality base calls.
                                char Base = read.Bases[baseIndex];
                                if (Base == variant.ReferenceAllele[0]) this.ReferenceCounts[varIndex]++;
                                if (Base == variant.VariantAlleles[0][0]) this.VariantCounts[varIndex]++;
                            }
                        }
                        break;
                    case 'S':
                        baseIndex += (int)cigar.Length;
                        break;
                    case 'I':
                        baseIndex += (int)cigar.Length;
                        break;
                    case 'D':
                        position += (int)cigar.Length;
                        break;
                    default:
                        // We don't know how to cope with this CIGAR operation; bail out!
                        return;
                }
            }
        }

        /// <summary>
        /// Step 3: Summarize results to a simple tab-delimited file and a CSV file.
        /// </summary>
        protected void WriteResults(string outputPath)
        {
            WriteAlleleCounts(outputPath);
            WriteBAlleleFrequencies(outputPath + ".baf");
        }

        protected void WriteAlleleCounts(string outputPath)
        {
            using (GzipWriter writer = new GzipWriter(outputPath))
            {
                writer.WriteLine("#Chromosome\tPosition\tRef\tAlt\tCountRef\tCountAlt");
                for (int index = 0; index < this.Variants.Count; index++)
                {
                    VcfVariant variant = this.Variants[index];
                    // skip HOM REF positions 
                    if (this.VariantCounts[index] > 5)
                    {
                        writer.WriteLine(string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", variant.ReferenceName, variant.ReferencePosition,
                        variant.ReferenceAllele, variant.VariantAlleles[0], this.ReferenceCounts[index],
                        this.VariantCounts[index]));
                    }
                }
            }
            Console.WriteLine("{0} Results written to {1}", DateTime.Now, outputPath);
        }

        protected void WriteBAlleleFrequencies(string outputPath)
        {
            using (FileStream stream = new FileStream(outputPath, FileMode.Create, FileAccess.Write))
            using (StreamWriter writer = new StreamWriter(stream))
            {
                writer.WriteLine(CSVWriter.GetLine("Chromosome", "Position", "BAF"));
                for (int index = 0; index < this.Variants.Count; index++)
                {
                    VcfVariant variant = this.Variants[index];
                    double? baf = GetBAlleleFrequency(variant, this.ReferenceCounts[index], this.VariantCounts[index]);
                    if (!baf.HasValue)
                        continue;
                    writer.WriteLine(CSVWriter.GetLine(variant.ReferenceName, variant.ReferencePosition.ToString(),
                        baf.Value.ToString()));
                }
            }
            Console.WriteLine("{0} Results written to {1}", DateTime.Now, outputPath);
        }

        public static double? GetBAlleleFrequency(VcfVariant variant, int referenceCount, int variantCount)
        {
            double? baf = null;
            double totalAlleleCount = referenceCount + variantCount;
            if (totalAlleleCount < 1)
                return baf;
            if (variant.ReferenceAllele.Equals(".") || variant.VariantAlleles[0].Equals("."))
                return baf;

            if (BAllelePreference(variant.ReferenceAllele) < BAllelePreference(variant.VariantAlleles[0]))
            {
                baf = referenceCount / totalAlleleCount;
            }
            else
            {
                baf = variantCount / totalAlleleCount;
            }

            return baf;
        }

        /// <summary>
        /// Returns B allele preference for single nucleotide alleles. The highest is 0 and the lowest is 3.
        /// This definition of B-Allele Frequency is based on the definition that is used for bead arrays.
        /// Here, the choice of the B allele is based on the color of dye attached to each nucleotide. A and T
        /// get one color, G and C get the other color. Bead array has much more complex rule for tie-breaking
        /// between A and T or G and C that involves top and bottom strands. This is unnecessary so we go with
        /// the simpler hierarchical approach. In this case, this is not strictly a definition of B allele
        /// frequency that is applied in any other product.
        /// </summary>
        /// <param name="allele"></param>
        /// <returns></returns>
        private static int BAllelePreference(string allele)
        {
            switch (allele.ToLower())
            {
                case "a":
                    return 0;
                case "t":
                    return 1;
                case "g":
                    return 2;
                case "c":
                    return 3;
                default:
                    throw new ArgumentException("Invalid single nucleotide allele: " + allele);
            }
        }
    }
}
