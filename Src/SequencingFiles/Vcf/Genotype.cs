using System;
using System.Collections.Generic;
using System.Linq;
using ILMNcommon.Common;

namespace SequencingFiles.Vcf
{
    public class Genotype
    {
        public List<Allele> Alleles { get; }

        public Genotype(List<Allele> alleles)
        {
            if (alleles.Empty())
                throw new ArgumentException("genotype cannot be empty");
            var nonNullAlleles = alleles.Where(allele => allele != null);
            ReferencePosition position = nonNullAlleles.Select(allele => allele.ReferenceStartPosition).FirstOrDefault();
            string reference = nonNullAlleles.Select(allele => allele.ReferenceAlleleSequence).FirstOrDefault();
            foreach (var allele in nonNullAlleles)
            {
                if (!position.Equals(allele.ReferenceStartPosition))
                    throw new ArgumentException($"All alleles must have the same reference start position. Found {position} and {allele.ReferenceStartPosition}");
                if (allele.ReferenceAlleleSequence != reference)
                    throw new ArgumentException($"All alleles must have the same reference allele sequence. Found reference sequences {reference} and {allele.ReferenceAlleleSequence} at {position}");
            }

            Alleles = alleles;
        }

        public static List<Genotype> GetGenotypes(VcfVariant variant)
        {
            var refStartPosition = new ReferencePosition(variant.ReferenceName, variant.ReferencePosition);
            Allele referenceAllele = Allele.CreateReference(refStartPosition, variant.ReferenceAllele);
            var alternateAlleles = Allele.GetAlternateAlleles(referenceAllele, variant.VariantAlleles);
            var genotypes = new List<Genotype>();
            foreach (var alleleIndexes in variant.GenotypeAlleleIndexes)
            {
                genotypes.Add(GetGenotype(referenceAllele, alternateAlleles, alleleIndexes));
            }
            return genotypes;
        }

        /// <summary>
        /// Create a genotype from the allele indexes specified in the GT field of a vcf entry
        /// </summary>
        /// <param name="referenceAllele">The reference allele</param>
        /// <param name="alternateAlleles">The alternate alleles from the ALT column of the Vcf entry</param>
        /// <param name="alleleIndexes">Allele indexes from the vcf GT field. Null means a missing allele (i.e. "."), 0 is the reference allele and allele index > 1 are the alternate alleles.</param>
        /// <returns></returns>
        public static Genotype GetGenotype(Allele referenceAllele, List<Allele> alternateAlleles, List<int?> alleleIndexes)
        {
            if (alleleIndexes == null) return null;
            List<Allele> genotype = new List<Allele>();
            foreach (var alleleIndex in alleleIndexes)
            {
                genotype.Add(GetHaplotype(referenceAllele, alternateAlleles, alleleIndex));
            }
            return new Genotype(genotype);
        }

        private static Allele GetHaplotype(Allele referenceAllele, List<Allele> altAlleles, int? alleleIndex)
        {
            if (alleleIndex == null) return null;
            if (alleleIndex == 0) return referenceAllele;
            if (alleleIndex < 0 || alleleIndex > altAlleles.Count)
                throw new ArgumentException($"Allele {alleleIndex.Value} is invalid. Must be positive and less than or equal to the number of ALT alleles {altAlleles.Count}.");

            return altAlleles[alleleIndex.Value - 1];
        }

        public bool IsVariant => Alleles.Where(haplotype => haplotype != null)
            .Any(haplotype => haplotype.AlleleType != AlleleType.Reference);

        /// <summary>
        /// True if this genotype has any alleles with one of the types in <see cref="alleleTypes"/>
        /// and only alleles of type <see cref="AlleleType.Reference"/> or one of the types in <see cref="alleleTypes"/>
        /// No alleles in this genotype may be null
        /// </summary>
        /// <param name="alleleTypes">The possible allele types</param>
        /// <returns></returns>
        public bool IsPure(params AlleleType[] alleleTypes)
        {
            List<AlleleType> alleleTypesAndRef = alleleTypes.ToList();
            alleleTypesAndRef.Add(AlleleType.Reference);
            return HasAny(alleleTypes) && HasOnly(alleleTypesAndRef.ToArray());
        }

        /// <summary>
        /// True if this genotype consists of only alleles with one of the types specified by <cref name="alleleTypes"/>
        /// No alleles in this genotype may be null
        /// </summary>
        /// <param name="alleleTypes">The possible allele types</param>
        /// <returns></returns>
        public bool HasOnly(params AlleleType[] alleleTypes)
        {
            return
                !IsMissing &&
                Alleles.TrueForAll(allele => alleleTypes.Contains(allele.AlleleType));
        }

        /// <summary>
        /// return true if this genotype contains all of the allele types specified by <cref name="alleleTypes"/>
        /// </summary>
        /// <param name="alleleTypes">The required allele types</param>
        /// <returns></returns>
        public bool HasAll(params AlleleType[] alleleTypes)
        {
            foreach (var alleleType in alleleTypes)
            {
                if (Alleles.All(allele => allele == null || allele.AlleleType != alleleType)) return false;
            }
            return true;
        }

        /// <summary>
        /// return true if this genotype contains any of the allele types specified by <cref name="alleleTypes"/>
        /// </summary>
        /// <param name="alleleTypes">The possible allele types</param>
        /// <returns></returns>
        public bool HasAny(params AlleleType[] alleleTypes)
        {
            if (Alleles.Any(allele => allele != null && alleleTypes.Contains(allele.AlleleType)))
                return true;
            return false;
        }

        /// <summary>
        /// True if all alleles are of type <see cref="AlleleType.Reference"/> 
        /// </summary>
        public bool IsReference => HasOnly(AlleleType.Reference);

        /// <summary>
        /// True if this genotype contains exactly two non-null alleles that have the same sequence
        /// </summary>
        public bool IsHomozygous
        {
            get
            {
                if (Alleles.Count != 2)
                    return false;
                if (Alleles[0] == null)
                    return false;
                return Alleles[0].Equals(Alleles[1]);
            }
        }

        /// <summary>
        /// True if any allele in this genotype is null
        /// </summary>
        public bool IsMissing => Alleles.Any(haplotype => haplotype == null);

        /// <summary>
        /// True if this genotype contains exactly two non-null alleles that have different sequences
        /// </summary>
        public bool IsHeterozygous
        {
            get
            {
                if (Alleles.Count != 2)
                    return false;
                if (IsMissing)
                    return false;
                return !Alleles[0].Equals(Alleles[1]);
            }
        }

        /// <summary>
        /// True if this genotype contains exactly one non-null allele
        /// </summary>
        public bool IsHemizygous
        {
            get
            {
                if (Alleles.Count != 1)
                    return false;
                if (IsMissing)
                    return false;
                return true;
            }
        }
    }
}