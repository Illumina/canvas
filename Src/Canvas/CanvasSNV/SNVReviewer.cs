using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using SequencingFiles;

namespace CanvasSNV
{
    /// <summary>
    /// The SNV reviewer takes as input a .vcf file (from the normal sample) and a .bam file (from the paired tumor sample).
    /// It processes a particular chromosome at a time.
    /// For all variants in the .vcf file that pass our filters, we review the variant allele frequency in the .bam file
    /// from the tumor sample.  We output that data in tabular format.
    /// </summary>
    class SNVReviewer
    {
        #region Members
        public string Chromosome;
        List<VcfVariant> Variants;
        int[] ReferenceCounts;
        int[] VariantCounts;
        int MinimumBaseQScore = 20; // You must be >= this base call quality to be counted.
        #endregion

        public int Main(string vcfPath, string bamPath, string outputPath)
        {
            
            if (!File.Exists(vcfPath))
            {
                Console.Error.WriteLine("Error: Input vcf file not found at {0}", vcfPath);
                return 1;
            }
            if (!File.Exists(bamPath))
            {
                Console.Error.WriteLine("Error: Input bam file not found at {0}", bamPath);
                return 1;
            }

            this.LoadVariants(vcfPath);
            ReferenceCounts = new int[Variants.Count];
            VariantCounts = new int[Variants.Count];
            this.ProcessBamFile(bamPath);
            this.WriteResults(outputPath);
            return 0;
        }

        /// <summary>
        /// Step 1: Load the normal het SNVs of interest.
        /// </summary>
        protected void LoadVariants(string vcfPath)
        {
            Console.WriteLine("{0} Loading variants of interest from {1}", DateTime.Now, vcfPath);
            this.Variants = new List<VcfVariant>();
            int overallCount = 0; 
            int countThisChromosome = 0;
            using (VcfReader reader = new VcfReader(vcfPath, requireGenotypes: false))
            {
                VcfVariant variant = new VcfVariant();
                while (true)
                {
                    bool result = reader.GetNextVariant(variant);
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
                    if ((variant.Genotypes != null && variant.Genotypes.Any()) && variant.Filters != "PASS") continue; // FILTER may not say PASS for a dbSNP VCF file
                    if (variant.Genotypes != null && variant.Genotypes.Any()) // not available if we use a dbSNP VCF file
                    {
                        if (!variant.Genotypes[0].ContainsKey("GT")) continue; // no genotype - we don't know if it's a het SNV.
                        string genotype = variant.Genotypes[0]["GT"];
                        if (genotype != "0/1" && genotype != "1/0") continue;

                        // Also require they have a high enough quality score:
                        if (variant.Genotypes[0].ContainsKey("GQX")) // Note: Allow no GQX field, in case we want to use another caller (e.g. Pisces) and not crash
                        {
                            float GQX = float.Parse(variant.Genotypes[0]["GQX"]);
                            if (GQX < 30) continue;
                        }
                    }
                    // Note: Let's NOT require the variant be in dbSNP.  Maybe we didn't do annotation, either because
                    // we chose not to or because we're on a reference without annotation available.
                    //if (variant.Identifier == ".") continue;
                    // Remember all the variants that pass all our tests:
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
                    if (read.RefID < 0 || read.RefID > refID) break; // We're past our chromosome of interest.
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
                    if (read.MapQuality == 0) continue;

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
                        for (int index = 0; index < cigar.Length; index++,position++,baseIndex++)
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
        /// Step 3: Summarize results to a simple tab-delimited file.
        /// </summary>
        protected void WriteResults(string outputPath)
        {
            using (GzipWriter writer = new GzipWriter(outputPath))
            {
                writer.WriteLine("#Chromosome\tPosition\tRef\tAlt\tCountRef\tCountAlt");
                for (int index = 0; index < this.Variants.Count; index++)
                { 
                    VcfVariant variant = this.Variants[index];
                    // skip HOM REF positions 
                    if (this.VariantCounts[index] > 5) {
                        writer.WriteLine(string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", variant.ReferenceName, variant.ReferencePosition,
                        variant.ReferenceAllele, variant.VariantAlleles[0], this.ReferenceCounts[index], 
                        this.VariantCounts[index]));
                    }
                }
            }
            Console.WriteLine("{0} Results written to {1}", DateTime.Now, outputPath);
        }

    }
}
