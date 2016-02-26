using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace SequencingFiles.Vcf
{
    /// <summary>
    ///     VcfVariant stores all of the data from one line of VCF data
    ///     The info fields and genotype columns are stored as strings
    /// </summary>
    public class VcfVariant
    {
        #region members
        public string Filters;
        public string[] GenotypeTagOrder = null;
        public List<Dictionary<string, string>> GenotypeColumns;
        public string Identifier;
        public Dictionary<string, string> InfoFields;
        public string[] InfoTagOrder = null;
        public double Quality;
        public bool HasQuality = true;
        public string ReferenceAllele;
        public string ReferenceName;
        public int ReferencePosition;
        public string[] VariantAlleles;
        public List<List<int?>> GenotypeAlleleIndexes;
        #endregion

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
            if (GenotypeColumns.Count <= sampleNumber) return false;
            var sampleGT = GenotypeColumns[sampleNumber];
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

            if (GenotypeColumns != null && GenotypeColumns.Any())
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
                foreach (Dictionary<string, string> t in GenotypeColumns)
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