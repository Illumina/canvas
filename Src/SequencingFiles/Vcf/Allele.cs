using System;
using System.Collections.Generic;
using System.Linq;
using ILMNcommon.Common;

namespace SequencingFiles.Vcf
{
    public enum AlleleType
    {
        Complex = 0, // A variant type we don't otherwise categorize - e.g. ref ATG, alt GTACC
        Reference = 1,
        Snv = 2,
        Insertion = 3,
        Deletion = 4,
        Mnv = 5, // multiple nucleotide variant
    }

    public class Allele : IEquatable<Allele>, IOverlappable
    {

        private readonly IOverlappable _overlapper;
        /// <summary>
        /// The index for this allele which would appear in the GT field of a VCF entry.
        /// 0 --> reference allele
        /// 1+ --> the one-based index of the alternate allele in the list of alternate alleles from the ALT column of a VCF entry 
        /// </summary>
        public int AlleleIndex { get; }
        public AlleleType AlleleType { get; set; }
        public ReferencePosition ReferenceStartPosition { get; }
        public string ReferenceAlleleSequence { get; }
        public string AlleleSequence { get; }

        /// <summary>
        /// <see cref="ReferenceAlleleSequence"/> after removing any common prefix or suffix with <see cref="AlleleSequence"/>.
        /// For a canonical insertion allele this will remove the padding base and return <see cref="string.Empty"/>
        /// </summary>
        public string TrimmedReferenceAlleleSequence { get; }

        /// <summary>
        /// <see cref="AlleleSequence"/> after removing any common prefix or suffix with <see cref="ReferenceAlleleSequence"/>.
        /// For a canonical deletion allele this will remove the padding base and return <see cref="string.Empty"/>
        /// </summary>
        public string TrimmedAlleleSequence { get; }

        private Allele(
            ReferencePosition referenceStartPosition,
            string referenceAlleleSequence,
            string alleleSequence,
            string trimmedReferenceAlleleSequence,
            string trimmedAlleleSequence,
            int alleleIndex,
            IOverlappable overlapper,
            AlleleType alleleType)
        {
            _overlapper = overlapper;
            if (string.IsNullOrWhiteSpace(referenceAlleleSequence))
                throw new ArgumentException("reference allele sequence cannot be null or empty");
            ReferenceStartPosition = referenceStartPosition;
            ReferenceAlleleSequence = referenceAlleleSequence;
            AlleleSequence = alleleSequence;
            TrimmedReferenceAlleleSequence = trimmedReferenceAlleleSequence;
            TrimmedAlleleSequence = trimmedAlleleSequence;
            AlleleIndex = alleleIndex;
            AlleleType = alleleType;
        }

        public static Allele CreateReference(ReferencePosition referenceStartPosition, string referenceAlleleSequence)
        {
            referenceAlleleSequence = referenceAlleleSequence.ToUpperInvariant();
            return new Allele(referenceStartPosition, referenceAlleleSequence, referenceAlleleSequence,
                "", "", 0, new NonVariantIntervalOverlapper(), AlleleType.Reference);
        }


        /// <summary>
        /// Create a variant allele from an alternate allele and the corresponding reference allele.
        /// The rules are as follows:
        /// - Trim off any common prefix or suffix.  Let |ref| denote the length of the
        ///   reference allele after trimming, and |alt| denote the length of the alt allele after trimming.
        /// - If |ref|=0, it's an insertion
        /// - If |alt|=0, it's a deletion
        /// - If |ref|=|alt|=1, it's a SNV
        /// - If |ref|=|alt| and |alt|>1, it's a MNP
        /// - Otherwise (i.e. |ref|>0 and |alt|>0 and |ref| != |alt|) it's a complex event
        /// </summary>
        /// <param name="referenceStartPosition"><see cref="ReferencePosition"/> corresponding to the first base in <see cref="referenceAllele"/></param>
        /// <param name="referenceAllele">The base sequence of the reference allele</param>
        /// <param name="alternateAllele">The base sequence of the alternate allele</param>
        /// <param name="alleleIndex">The one-based index of <see cref="alternateAllele"/> in the list of alternate alleles from the ALT column of the VCF entry</param>
        /// <returns></returns>
        public static Allele CreateVariant(ReferencePosition referenceStartPosition, string referenceAllele, string alternateAllele, int alleleIndex)
        {
            if (alleleIndex <= 0)
                throw new ArgumentException("allele index for a variant allele must be > 0");
            if (referenceAllele == null)
                throw new ArgumentException("reference allele cannot be null");
            if (alternateAllele == null)
                throw new ArgumentException("alternate allele cannot be null");
            referenceAllele = referenceAllele.ToUpperInvariant();
            alternateAllele = alternateAllele.ToUpperInvariant();
            if (referenceAllele == alternateAllele)
                throw new ArgumentException("alternate allele cannot be the same as reference allele");
            if (alternateAllele == ".")
                throw new ArgumentException("alternate allele cannot be \".\"");

            int commonSuffixBases = referenceAllele.CommonSuffixLength(alternateAllele);
            string trimmedReferenceAllele = referenceAllele.Substring(0, referenceAllele.Length - commonSuffixBases);
            string trimmedAlternateAllele = alternateAllele.Substring(0, alternateAllele.Length - commonSuffixBases);
            int commonPrefixBases = trimmedReferenceAllele.CommonPrefixLength(trimmedAlternateAllele);
            trimmedReferenceAllele = trimmedReferenceAllele.Substring(commonPrefixBases);
            trimmedAlternateAllele = trimmedAlternateAllele.Substring(commonPrefixBases);

            ReferencePosition refStartPositionAfterTrim = referenceStartPosition.Shift(commonPrefixBases);

            IOverlappable overlapper;
            AlleleType alleleType;
            if (trimmedReferenceAllele.Length == 0)
            {
                alleleType = AlleleType.Insertion;
                overlapper = new InsertionOverlapper(refStartPositionAfterTrim.Previous(), referenceAllele.Substring(commonPrefixBases), trimmedAlternateAllele);
            }
            else
            {
                // this variant could potentially start at multiple positions along the reference (e.g. ref TT and alt T) and still result in the same sequence
                // we want to allow for this when checking for overlapping variant bases
                ReferencePosition rightmostRefStartPositionAfterTrim = referenceStartPosition.Shift(referenceAllele.CommonPrefixLength(alternateAllele));
                ReferenceInterval referenceIntervalAfterTrim = new ReferenceInterval(refStartPositionAfterTrim, rightmostRefStartPositionAfterTrim.Shift(trimmedReferenceAllele.Length - 1));
                overlapper = new VariantIntervalOverlapper(referenceIntervalAfterTrim);
                if (trimmedAlternateAllele.Length == 0)
                    alleleType = AlleleType.Deletion;
                else if (trimmedAlternateAllele.Length == 1 && trimmedReferenceAllele.Length == 1)
                    alleleType = AlleleType.Snv;
                else if (trimmedAlternateAllele.Length == trimmedReferenceAllele.Length)
                    alleleType = AlleleType.Mnv;
                else
                    alleleType = AlleleType.Complex;
            }

            return new Allele(referenceStartPosition, referenceAllele, alternateAllele,
                trimmedReferenceAllele, trimmedAlternateAllele, alleleIndex, overlapper, alleleType);
        }

        public bool VariantBasesOverlap(ReferenceInterval interval)
        {
            return _overlapper.VariantBasesOverlap(interval);
        }

        public static List<Allele> GetAlternateAlleles(VcfVariant vcfEntry)
        {
            var referenceAllele = GetReferenceAllele(vcfEntry);
            return GetAlternateAlleles(referenceAllele, vcfEntry.VariantAlleles);
        }

        public static Allele GetReferenceAllele(VcfVariant vcfEntry)
        {
            ReferencePosition referenceStartPosition = new ReferencePosition(vcfEntry.ReferenceName, vcfEntry.ReferencePosition);
            return CreateReference(referenceStartPosition, vcfEntry.ReferenceAllele);
        }

        public static List<Allele> GetAlternateAlleles(Allele referenceAllele, IEnumerable<string> alternateAlleleSequences)
        {
            var alternateAlleles = new List<Allele>();
            if (alternateAlleleSequences.Count() == 1 && alternateAlleleSequences.First() == ".")
                return alternateAlleles;

            int alleleIndex = 1;
            foreach (var altHaplotypeSequence in alternateAlleleSequences)
            {
                var altHaplotype = CreateVariant(referenceAllele.ReferenceStartPosition, referenceAllele.ReferenceAlleleSequence, altHaplotypeSequence, alleleIndex);
                alternateAlleles.Add(altHaplotype);
                alleleIndex++;
            }
            return alternateAlleles;
        }

        public bool IsTransition
        {
            get
            {
                if (AlleleType != AlleleType.Snv) return false;
                char referenceBase = TrimmedReferenceAlleleSequence.Single();
                char variantBase = TrimmedAlleleSequence.Single();
                // no need to convert to upper here since the static factory methods ensure uppercase
                return (((referenceBase == 'A') && (variantBase == 'G')) ||
                        ((referenceBase == 'G') && (variantBase == 'A')) ||
                        ((referenceBase == 'C') && (variantBase == 'T')) ||
                        ((referenceBase == 'T') && (variantBase == 'C')));
            }
        }

        public bool IsTransversion
        {
            get
            {
                if (AlleleType != AlleleType.Snv) return false;
                return !IsTransition;
            }
        }

        #region IEquatable

        public override bool Equals(object obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;

            Allele allele = obj as Allele;
            if (allele == null)
                return false;

            return Equals(allele);
        }

        public bool Equals(Allele allele)
        {
            if (ReferenceEquals(null, allele)) return false;
            if (ReferenceEquals(this, allele)) return true;

            return
                ReferenceStartPosition.Equals(allele.ReferenceStartPosition) &&
                ReferenceAlleleSequence == allele.ReferenceAlleleSequence &&
                AlleleSequence == allele.AlleleSequence;
        }

        public override int GetHashCode()
        {
            int hash = 23;
            hash = hash * 31 + ReferenceStartPosition.GetHashCode();
            hash = hash * 31 + ReferenceAlleleSequence.GetHashCode();
            return hash * 31 + AlleleSequence.GetHashCode();
        }

        #endregion
    }
}