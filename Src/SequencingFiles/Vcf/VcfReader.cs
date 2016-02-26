using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using SequencingFiles.Vcf;

namespace SequencingFiles
{
    /// <summary>
    ///     VcfReader is a multi-sample vcf reader.
    /// </summary>
    public class VcfReader : IDisposable
    {
        #region members
        public List<string> HeaderLines = new List<string>();
        private bool IsDisposed;
        private bool IsOpen;
        private GzipReader Reader;
        public List<string> Samples = new List<string>();
        private static char[] InfoSplitChars = new char[] { ';' };
        private bool RequireGenotypes;
        // For (minor) speedup, cache the genotype tag order from one line to the next, because it's typically the same for every record:
        private string GenotypeTagString;
        private string[] GenotypeTagOrder;
        #endregion

        // constructor
        public VcfReader(string vcfPath, bool requireGenotypes = true, bool skipHeader = false)
        {
            this.RequireGenotypes = requireGenotypes;
            IsOpen = false;
            Open(vcfPath, skipHeader);
        }

        #region IDisposable
        // Note: These two pages explain IDisposable in great detail and give a picture for Why We Do Things This Way:
        // http://stackoverflow.com/questions/538060/proper-use-of-the-idisposable-interface
        // http://msdn.microsoft.com/en-us/library/system.idisposable(v=vs.110).aspx
        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        // destructor
        ~VcfReader()
        {
            Dispose(false);
        }

        protected virtual void Dispose(bool disposing)
        {
            lock (this)
            {
                if (!IsDisposed)
                {
                    IsDisposed = true;
                    Close();
                }
            }
        }
        #endregion

        /// <summary>
        ///     populates a vcf variant object given an array of vcf columns
        /// </summary>
        protected void ConvertColumnsToVariant(string[] cols, VcfVariant variant)
        {
            variant.ReferenceName = cols[VcfCommon.ChromIndex];
            variant.ReferencePosition = int.Parse(cols[VcfCommon.PosIndex]);
            variant.Identifier = cols[VcfCommon.IDIndex];
            variant.ReferenceAllele = cols[VcfCommon.RefIndex];
            variant.Filters = cols[VcfCommon.FilterIndex];

            if (cols[VcfCommon.QualIndex] == ".")
                variant.HasQuality = false;
            double.TryParse(cols[VcfCommon.QualIndex], out variant.Quality); // CFTR uses a ".", which is not actually legal... (actually, vcf 4.1 does allow the missing value "." here. Strelka uses it)

            // parse the variant alleles
            variant.VariantAlleles = cols[VcfCommon.AltIndex].Split(',');

            // parse the info fields
            //variant.InfoFields.Clear();
            variant.InfoFields = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
            string InfoData = cols[VcfCommon.InfoIndex];
            if (InfoData == ".") InfoData = ""; // Special case: a "." in the INFO field should be treated like an empty string.
            string[] infoCols = InfoData.Split(InfoSplitChars, StringSplitOptions.RemoveEmptyEntries);

            int numInfoCols = infoCols.Length;

            if ((variant.InfoTagOrder == null) || (numInfoCols != variant.InfoTagOrder.Length))
            {
                variant.InfoTagOrder = new string[numInfoCols];
            }

            for (int infoColIndex = 0; infoColIndex < numInfoCols; infoColIndex++)
            {
                string infoField = infoCols[infoColIndex];
                string[] infoFieldKvp = infoField.Split('=');
                variant.InfoTagOrder[infoColIndex] = infoFieldKvp[0];
                variant.InfoFields[infoFieldKvp[0]] = (infoFieldKvp.Length == 1 ? null : infoFieldKvp[1]);
            }

            if (cols.Length > VcfCommon.GenotypeIndex) // Genotype columns present
            {
                // parse the genotype format field
                if (cols[VcfCommon.FormatIndex] != GenotypeTagString)
                {
                    GenotypeTagString = cols[VcfCommon.FormatIndex];
                    GenotypeTagOrder = GenotypeTagString.Split(':');
                }
                variant.GenotypeTagOrder = GenotypeTagOrder;

                // parse the genotype data for each sample
                variant.GenotypeColumns = new List<Dictionary<string, string>>();
                for (int sampleIndex = 0; sampleIndex < this.Samples.Count; sampleIndex++)
                {
                    string genotypeColumn = cols[VcfCommon.GenotypeIndex + sampleIndex];
                    if (genotypeColumn == ".")
                    {
                        variant.GenotypeColumns.Add(null);
                    }
                    else
                    {
                        string[] genotypeCols = genotypeColumn.Split(':');
                        variant.GenotypeColumns.Add(ParseGenotype(variant.GenotypeTagOrder, genotypeCols));
                    }
                }
                variant.GenotypeAlleleIndexes = GetGenotypeAlleleIndexes(variant);
            }
        }

        private static List<List<int?>> GetGenotypeAlleleIndexes(VcfVariant variant)
        {
            var sampleAlleleIndexes = new List<List<int?>>();
            foreach (var genotype in variant.GenotypeColumns)
            {
                if (genotype == null || !genotype.ContainsKey(VcfCommon.GenotypeColumnFields.Genotype))
                    sampleAlleleIndexes.Add(null);
                else
                    sampleAlleleIndexes.Add(GetAlleleIndexes(genotype[VcfCommon.GenotypeColumnFields.Genotype]));
            }
            return sampleAlleleIndexes;
        }

        private static List<int?> GetAlleleIndexes(string genotype)
        {
            List<int?> alleleIndexes = new List<int?>();
            foreach (var unparsedHaplotype in genotype.Split('/', '|'))
            {
                alleleIndexes.Add(ParseAlleleIndex(unparsedHaplotype));
            }
            return alleleIndexes;
        }

        private static int? ParseAlleleIndex(string unparsedHaplotype)
        {
            if (unparsedHaplotype == ".")
                return null;

            uint alleleIndex;
            if (!uint.TryParse(unparsedHaplotype, out alleleIndex))
                throw new FormatException($"{VcfCommon.GenotypeColumnFields.Genotype} field contains invalid string {unparsedHaplotype}.");
            return (int)alleleIndex;
        }

        /// <summary>
        ///     closes the vcf file
        /// </summary>
        private void Close()
        {
            if (!IsOpen) return;
            IsOpen = false;
            Reader.Close();
        }

        /// <summary>
        /// Loop over variants like this: foreach (VcfVariant variant in reader.GetVariants())
        /// </summary>
        public IEnumerable<VcfVariant> GetVariants()
        {
            // sanity check: make sure the file is open
            if (!IsOpen) yield break;

            while (true)
            {
                // grab the next vcf line
                string line = Reader.ReadLine();
                if (line == null) break;

                VcfVariant variant = new VcfVariant();

                // split the columns and assign them to VcfVariant
                string[] cols = line.Split('\t');

                // convert the columns to a variant
                ConvertColumnsToVariant(cols, variant);
                if (RequireGenotypes && (variant.GenotypeColumns == null || variant.GenotypeColumns.Count == 0))
                    throw new ApplicationException("Missing genotype columns in VCF file");
                yield return variant;
            }
        }

        /// <summary>
        /// Test method: Load variant but keep unparsed line around
        /// </summary>
        public bool GetNextVariant(VcfVariant variant, out string line)
        {
            line = null;
            // sanity check: make sure the file is open
            if (!IsOpen) return false;

            // grab the next vcf line
            line = Reader.ReadLine();
            if (line == null) return false;

            // split the columns and assign them to VcfVariant
            string[] cols = line.Split('\t');

            // convert the columns to a variant
            ConvertColumnsToVariant(cols, variant);
            if (RequireGenotypes && variant.GenotypeColumns.Count == 0)
                throw new ApplicationException("Missing genotype columns in VCF file");
            return true;
        }

        /// <summary>
        /// Test method: Load variant but keep unparsed line around
        /// </summary>
        public bool GetNextVariant(VcfVariant variant, out string line, out string[] bits)
        {
            line = null;
            bits = null;
            // sanity check: make sure the file is open
            if (!IsOpen) return false;

            // grab the next vcf line
            line = Reader.ReadLine();
            if (line == null) return false;

            // split the columns and assign them to VcfVariant
            bits = line.Split('\t');

            // convert the columns to a variant
            ConvertColumnsToVariant(bits, variant);
            if (RequireGenotypes && variant.GenotypeColumns.Count == 0)
                throw new ApplicationException("Missing genotype columns in VCF file");
            return true;
        }

        /// <summary>
        ///     Retrieves the next available variant and returns false if no variants are available.
        /// </summary>
        public bool GetNextVariant(VcfVariant variant)
        {
            // sanity check: make sure the file is open
            if (!IsOpen) return false;

            // grab the next vcf line
            string line = Reader.ReadLine();
            if (line == null) return false;

            // split the columns and assign them to VcfVariant
            string[] cols = line.Split('\t');

            // convert the columns to a variant
            ConvertColumnsToVariant(cols, variant);
            if (RequireGenotypes && variant.GenotypeColumns.Count == 0)
                throw new ApplicationException("Missing genotype columns in VCF file");

            return true;
        }

        /// <summary>
        ///     opens the vcf file and reads the header
        /// </summary>
        private void Open(string vcfPath, bool skipHeader)
        {
            // sanity check: make sure the vcf file exists
            if (!File.Exists(vcfPath))
            {
                throw new FileNotFoundException(string.Format("The specified vcf file ({0}) does not exist.", vcfPath));
            }

            Reader = new GzipReader(vcfPath);
            IsOpen = true;
            if (skipHeader)
            {
                this.Samples.Add("Sample");
            }
            else
            {
                ParseHeader();
            }
        }

        /// <summary>
        ///     parse a sample genotype column and returns the corresponding dictionary
        /// </summary>
        private static Dictionary<string, string> ParseGenotype(string[] genotypeFormatTags, string[] genotypeCols)
        {
            Dictionary<string, string> genotypeMap = new Dictionary<string, string>();
            // sanity check: make sure we have the same number of columns
            if (genotypeFormatTags.Length < genotypeCols.Length)
            {
                throw new ApplicationException(string.Format(
                        "VCF parse error: Expected the same number of columns in the genotype format column ({0}) as in the sample genotype column ({1}).",
                        genotypeFormatTags.Length, genotypeCols.Length));
            }
            for (int colIndex = 0; colIndex < genotypeCols.Length; colIndex++)
            {
                genotypeMap[genotypeFormatTags[colIndex]] = genotypeCols[colIndex].Trim().Replace("\"", "");
            }

            return genotypeMap;
        }

        /// <summary>
        ///     reads the vcf header
        /// </summary>
        private void ParseHeader()
        {
            // store the header
            string line;

            while (true)
            {
                // grab the next line - stop if we have reached the main header or read the entire file
                line = Reader.ReadLine();
                if ((line == null) || line.StartsWith(VcfCommon.ChromosomeHeader)) break;
                HeaderLines.Add(line);
            }

            // sanity check
            if ((line == null) || !line.StartsWith(VcfCommon.ChromosomeHeader))
            {
                throw new ApplicationException(
                    string.Format("Could not find the vcf header (starts with {0}). Is this a valid vcf file?",
                                  VcfCommon.ChromosomeHeader));
            }

            // establish how many samples we have
            string[] headerCols = line.Split('\t');
            HeaderLines.Add(line);
            int sampleCount = headerCols.Length - VcfCommon.GenotypeIndex;

            for (int sampleIndex = 0; sampleIndex < sampleCount; sampleIndex++)
            {
                Samples.Add(headerCols[VcfCommon.GenotypeIndex + sampleIndex]);
            }
        }

        /// <summary>
        /// Load a list of all variants in a file.  This is memory-intensive; don't do this for whole-genome vcf files!
        /// </summary>
        public static List<VcfVariant> GetAllVariantsInFile(string vcfPath)
        {
            using (VcfReader reader = new VcfReader(vcfPath))
            {
                return reader.GetVariants().ToList();
            }
        }

        /// <summary>
        ///     Returns the actual position within the vcf file (used when making our ad-hoc index)
        /// </summary>
        public long Position()
        {
            return Reader.GetCurrentPosition();
        }
    }
}