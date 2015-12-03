using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace SequencingFiles
{
    /// <summary>
    ///     VcfVariant stores all of the data from one line of VCF data
    ///     Let's keep the info fields and genotypes stored as strings
    ///     for now.
    /// </summary>
    public class VcfVariant
    {
        #region members
        public string Filters;
        public string[] GenotypeTagOrder = null;
        public List<Dictionary<string, string>> Genotypes;
        public string Identifier;
        public Dictionary<string, string> InfoFields;
        public string[] InfoTagOrder = null;
        public double Quality;
        public bool HasQuality = true;
        public string ReferenceAllele;
        public string ReferenceName;
        public int ReferencePosition;
        public VariantType VarType1; // Variant type of allele 1
        public VariantType VarType2; // Variant type of allele 2
        public string[] VariantAlleles;
        #endregion

        public static VariantType GetVariantType(string referenceAllele, string variantAllele)
        {
            if (referenceAllele == variantAllele) return VariantType.Reference;
            if (referenceAllele.Length == 1 && variantAllele.Length == 1) return VariantType.SNV;
            if (referenceAllele.Length > variantAllele.Length && referenceAllele.StartsWith(variantAllele))
                return VariantType.Deletion;
            if (referenceAllele.Length < variantAllele.Length && variantAllele.StartsWith(referenceAllele))
                return VariantType.Insertion;
            return VariantType.Complex; // Complex!
        }

        /// <summary>
        /// Return a description of this .vcf record type
        /// </summary>
        public VariantType GetOverallVarTypeDeprecated()
        {
            switch (this.VarType1)
            {
                case VariantType.Reference:
                    return this.VarType2;
                case VariantType.SNV:
                    switch (this.VarType2)
                    {
                        case VariantType.Reference:
                        case VariantType.SNV:
                            return VariantType.SNV;
                        case VariantType.Insertion:
                            return VariantType.SNVInsertion;
                        case VariantType.Deletion:
                            return VariantType.SNVDeletion;
                        default:
                            return VariantType.Complex;
                    }
                case VariantType.MNP:
                    switch (this.VarType2)
                    {
                        case VariantType.Reference:
                        case VariantType.MNP:
                            return VariantType.MNP;
                        default:
                            return VariantType.Complex;
                    }
                case VariantType.Insertion:
                    switch (this.VarType2)
                    {
                        case VariantType.Reference:
                        case VariantType.Insertion:
                            return VariantType.Insertion;
                        case VariantType.SNV:
                            return VariantType.SNVInsertion;
                        case VariantType.Deletion:
                            return VariantType.InsertionDeletion;
                        default:
                            return VariantType.Complex;
                    }
                case VariantType.Deletion:
                    switch (this.VarType2)
                    {
                        case VariantType.Reference:
                        case VariantType.Deletion:
                            return VariantType.Deletion;
                        case VariantType.SNV:
                            return VariantType.SNVDeletion;
                        case VariantType.Insertion:
                            return VariantType.InsertionDeletion;
                        default:
                            return VariantType.Complex;
                    }
                default:
                    return VariantType.Complex;
            }
        }


        /// <summary>
        /// Return TRUE if this vcf record is a homozygous or hemizygous reference call (not a variant, not a no-call).
        /// </summary>
        public bool IsReferenceCall()
        {
            //if (this.VarType == VariantType.Missing) return false;
            return this.VarType1 == VariantType.Reference && (this.VarType2 == VariantType.Missing || this.VarType2 == VariantType.Reference);
        }

        /// <summary>
        /// Add a filter, checking to make sure we don't add it twice:
        /// </summary>
        public void AddFilter(string filter)
        {
            if (this.Filters == "PASS")
            {
                this.Filters = filter;
                return;
            }
            if (this.Filters.Split(';').Contains(filter)) return;
            this.Filters = this.Filters + ";" + filter;
        }

        public bool IsSNV(int genotype1, int genotype2) // %%% DEPRECATED - use VarType
        {
            if (ReferenceAllele.Length != 1) return false;
            if (genotype1 > 0 && VariantAlleles[genotype1 - 1].Length == 1) return true;
            if (genotype2 > 0 && VariantAlleles[genotype2 - 1].Length == 1) return true;
            return false;
        }

        public bool IsDeletion(int genotype1, int genotype2) // %%% DEPRECATED - use VarType
        {
            if (genotype1 > 0 && VariantAlleles[genotype1 - 1].Length < ReferenceAllele.Length) return true;
            if (genotype2 > 0 && VariantAlleles[genotype2 - 1].Length < ReferenceAllele.Length) return true;
            return false;
        }

        public bool IsInsertion(int genotype1, int genotype2) // %%% DEPRECATED - use VarType
        {
            if (genotype1 > 0 && VariantAlleles[genotype1 - 1].Length > ReferenceAllele.Length) return true;
            if (genotype2 > 0 && VariantAlleles[genotype2 - 1].Length > ReferenceAllele.Length) return true;
            return false;
        }

        public int LengthVariant(int variantAlleleIndex)
        {
            int varTypeInt = (int) GetVariantType(ReferenceAllele, VariantAlleles[variantAlleleIndex]);
            string altAllele = VariantAlleles[variantAlleleIndex];

            switch (varTypeInt)
            {
                case 2: //SNP
                    return 1;

                case 4: //if IsDeletion
                    string dif1 = ReferenceAllele.Substring(altAllele.Length);
                    return dif1.Length;

                case 3: //if IsInsertion
                    string dif2 = altAllele.Substring(ReferenceAllele.Length); // get the substring 
                    return dif2.Length;

                case 0: //if IsComplexVariant, this is the best we can do for an unknown type of variant
                    return altAllele.Length;

                default:
                    throw new ApplicationException("The given genotype is not supported for this method"); // ???
            }
        }

        public string GetCall(int genotype1, int genotype2)
        {
            return GetCall(genotype1, genotype2, false);
        }

        /// <summary>
        ///     There is many functions use GetCall(int Genotype1, int Genotype2), so make a new one with
        ///     bool withBase for only variants display in the Variants table. See bug 59346
        /// </summary>
        /// <param name="genotype1"></param>
        /// <param name="genotype2"></param>
        /// <param name="withBase"></param>
        /// <returns></returns>
        public string GetCall(int genotype1, int genotype2, bool withBase)
        {
            if (IsSNV(genotype1, genotype2))
            {
                if (genotype1 == 0)
                {
                    return string.Format("{0}->{0}{1}", ReferenceAllele, VariantAlleles[genotype2 - 1]);
                }
                if (genotype2 == 0)
                {
                    return string.Format("{0}->{0}{1}", ReferenceAllele, VariantAlleles[genotype1 - 1]);
                }
                return string.Format("{0}->{1}{2}", ReferenceAllele, VariantAlleles[genotype1 - 1],
                    VariantAlleles[genotype2 - 1]);
            }
            StringBuilder sb = new StringBuilder();
            if (IsDeletion(genotype1, genotype2))
            {
                string altAllele = (genotype1 != 0 ? VariantAlleles[genotype1 - 1] : VariantAlleles[genotype2 - 1]);
                if (altAllele.Length >= ReferenceAllele.Length)
                {
                    altAllele = (genotype2 != 0 ? VariantAlleles[genotype2 - 1] : VariantAlleles[genotype1 - 1]);
                }
                string dif = ReferenceAllele.Substring(altAllele.Length);
                for (int i = 0; i < dif.Length; i++)
                    sb.Append("-");

                return ((withBase)
                    ? string.Format("{0}{1}/{0}{2}", altAllele, dif, sb)
                    : string.Format("{0}/{1}", dif, sb));
            }
            if (IsInsertion(genotype1, genotype2))
            {
                string altAllele = (genotype1 != 0 ? VariantAlleles[genotype1 - 1] : VariantAlleles[genotype2 - 1]);
                if (altAllele.Length <= ReferenceAllele.Length)
                {
                    altAllele = (genotype2 != 0 ? VariantAlleles[genotype2 - 1] : VariantAlleles[genotype1 - 1]);
                }

                string dif = altAllele.Substring(ReferenceAllele.Length); // get the substring 
                for (int i = 0; i < dif.Length; i++)
                    sb.Append("-");

                return ((withBase)
                    ? string.Format("{0}{1}/{0}{2}", ReferenceAllele, sb, dif)
                    : string.Format("{0}/{1}", sb, dif));
            }
            return ""; // ???
        }

        #region utility functions for the info field

        /// <summary>
        ///     returns true and sets the value if the key exists
        /// </summary>
        public bool TryParseInfoDouble(string key, out double value)
        {
            value = 0;
            string valueString;
            if (!InfoFields.TryGetValue(key, out valueString)) return false;
            if (!double.TryParse(valueString, out value)) return false;
            return true;
        }

        /// <summary>
        ///     returns true and sets the value if the key exists
        /// </summary>
        public bool TryParseInfoInt(string key, out int value)
        {
            value = 0;
            string valueString;
            if (!InfoFields.TryGetValue(key, out valueString)) return false;
            if (!int.TryParse(valueString, out value)) return false;
            return true;
        }

        /// <summary>
        ///     Parse values from sample genotype. Return true if the key exists.
        ///     sampleNumber is 0 based
        /// </summary>
        public bool TryGetGenotypeField(int sampleNumber, string key, out string value)
        {
            value = null;
            if (Genotypes.Count <= sampleNumber) return false;
            var sampleGT = Genotypes[sampleNumber];
            if (!sampleGT.TryGetValue(key, out value)) return false;
            return true;
        }

        #endregion

        /// <summary>
        /// Slightly awkward code to add (or overwrite) an INFO tag:
        /// </summary>
        public void AddInfoTag(string Key, string Value)
        {
            if (!this.InfoTagOrder.Contains(Key))
            {
                string[] newTagOrder = new string[this.InfoTagOrder.Length + 1];
                Array.Copy(this.InfoTagOrder, newTagOrder, this.InfoTagOrder.Length);
                newTagOrder[newTagOrder.Length - 1] = Key;
                this.InfoTagOrder = newTagOrder;
            }
            this.InfoFields[Key] = Value;
        }

        #region ToString

        /// <summary>
        ///     Re-encode the information in the vcf variant into a vcf line
        /// </summary>
        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendFormat("{0}\t{1}\t{2}\t{3}\t", ReferenceName, ReferencePosition, Identifier, ReferenceAllele);

            // write the alternate alleles
            int numVariantAlleles = VariantAlleles.Length;
            sb.Append(VariantAlleles[0]);

            for (int variantIndex = 1; variantIndex < numVariantAlleles; variantIndex++)
            {
                sb.AppendFormat(",{0}", VariantAlleles[variantIndex]);
            }

            if (HasQuality)
                sb.AppendFormat("\t{0:0.##}", Quality);
            else
                sb.Append("\t.");

            sb.AppendFormat("\t{0}\t", Filters);

            // write the info column
            bool needSemicolon = false;
            foreach (string infoField in InfoTagOrder)
            {
                if (needSemicolon) sb.Append(';');

                string infoValue;
                if (!InfoFields.TryGetValue(infoField, out infoValue))
                {
                    throw new ApplicationException(
                        string.Format("Unable to find the dictionary entry for the following info key: {0}", infoField));
                }

                if (infoValue == null)
                {
                    sb.Append(infoField);
                }
                else
                {
                    sb.AppendFormat("{0}={1}", infoField, infoValue);
                }

                needSemicolon = true;
            }
            if (InfoTagOrder.Length == 0) sb.Append("."); // VCF spec says: Any empty fields must be '.' 

            if (Genotypes != null && Genotypes.Any()) 
            {
                sb.Append('\t');

                // write the format field
                bool needColon = false;
                foreach (string formatField in GenotypeTagOrder)
                {
                    if (needColon) sb.Append(':');
                    sb.Append(formatField);
                    needColon = true;
                }

                // write the sample genotypes
                foreach (Dictionary<string, string> t in Genotypes)
                {
                    sb.Append('\t');
                    if (t == null)
                    {
                        sb.Append(".");
                    }
                    else
                    {
                        Dictionary<string, string> currentGenotype = t;

                        needColon = false;
                        foreach (string formatField in GenotypeTagOrder)
                        {
                            if (needColon) sb.Append(':');

                            string formatValue;
                            if (!currentGenotype.TryGetValue(formatField, out formatValue))
                            {
                                throw new ApplicationException(
                                    string.Format("Unable to find the dictionary entry for the following info key: {0}",
                                        formatField));
                            }

                            sb.Append(formatValue);
                            needColon = true;
                        }
                    }
                }
            }

            return sb.ToString();
        }

        #endregion


        #region Batch operation
        public delegate T VcfOperation<out T>(VcfVariant variant);
        public static List<T> OperateOnVariantsInFile<T>(string fileName, VcfOperation<T> operation)
        {
            if (!File.Exists(fileName))
                return null;
            List<T> variantList = new List<T>();

            using (VcfReader reader = new VcfReader(fileName))
            {
                foreach (VcfVariant variant in reader.GetVariants())
                {
                    variantList.Add(operation(variant));
                }
            }
            return variantList;
        }

        #endregion

    }
}