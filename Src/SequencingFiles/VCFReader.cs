using System;
using System.Collections.Generic;
using System.IO;
using System.Text;

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
			variant.InfoFields = new Dictionary<string, string>();
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
				string[] infoFieldKVP = infoField.Split('=');

				variant.InfoTagOrder[infoColIndex] = infoFieldKVP[0];
				variant.InfoFields[infoFieldKVP[0]] = (infoFieldKVP.Length == 1 ? null : infoFieldKVP[1]);
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
				variant.Genotypes = new List<Dictionary<string, string>>();
				for (int sampleIndex = 0; sampleIndex < this.Samples.Count; sampleIndex++)
				{
					string genotypeColumn = cols[VcfCommon.GenotypeIndex + sampleIndex];
					if (genotypeColumn == ".")
					{
						variant.Genotypes.Add(null);
					}
					else
					{
						string[] genotypeCols = genotypeColumn.Split(':');
						variant.Genotypes.Add(ParseGenotype(variant.GenotypeTagOrder, genotypeCols));
					}
				}

				// specify the variant type:
				AssignVariantType(variant);
			}
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
				if (RequireGenotypes && (variant.Genotypes == null || variant.Genotypes.Count == 0))
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
			if (RequireGenotypes && variant.Genotypes.Count == 0)
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
			if (RequireGenotypes && variant.Genotypes.Count == 0)
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
			if (RequireGenotypes && variant.Genotypes.Count == 0)
				throw new ApplicationException("Missing genotype columns in VCF file");

			return true;
		}

		private static void AssignVariantType(VcfVariant variant)
		{
			string genotype = null;

			if (variant.Genotypes[0] != null && variant.Genotypes[0].ContainsKey("GT"))
			{
				genotype = variant.Genotypes[0]["GT"];
			}

			// sanity check: support missing genotypes
			if (genotype == null || genotype == "./." || genotype == ".")
			{
				variant.VarType1 = VariantType.Missing;
				variant.VarType2 = VariantType.Missing;
				variant.VarType = VariantType.Missing;
				return;
			}
			// Handle usual cases like 0/0, 0/1, 1/0, 1/1 as well as 
			// special cases like ., ./., ./1, 1/.:
			int haplotypeA = int.TryParse(genotype.Substring(0, 1), out haplotypeA) ? haplotypeA : -1;
			int haplotypeB = genotype.Length >= 3 && int.TryParse(genotype.Substring(2, 1), out haplotypeB) ? haplotypeB : -1;
			// Treat things like ./1 or 0/. as homozygous:
			if (haplotypeA == -1) haplotypeA = haplotypeB;
			if (haplotypeB == -1) haplotypeB = haplotypeA;

			variant.VarType1 = GetAlleleVariantType(variant, haplotypeA);
			variant.VarType2 = GetAlleleVariantType(variant, haplotypeB);

			switch (variant.VarType1)
			{
				case VariantType.Reference:
					variant.VarType = variant.VarType2;
					break;
				case VariantType.SNV:
					switch (variant.VarType2)
					{
						case VariantType.Reference:
							variant.VarType = VariantType.SNV;
							break;
						case VariantType.SNV:
							variant.VarType = VariantType.SNV;
							break;
						case VariantType.Insertion:
							variant.VarType = VariantType.SNVInsertion;
							break;
						case VariantType.Deletion:
							variant.VarType = VariantType.SNVDeletion;
							break;
						default:
							variant.VarType = VariantType.Complex;
							break;
					}
					break;
				case VariantType.MNP:
					switch (variant.VarType2)
					{
						case VariantType.Reference:
							variant.VarType = VariantType.MNP;
							break;
						case VariantType.MNP:
							variant.VarType = VariantType.MNP;
							break;
						default:
							variant.VarType = VariantType.Complex;
							break;
					}
					break;
				case VariantType.Insertion:
					switch (variant.VarType2)
					{
						case VariantType.Reference:
							variant.VarType = VariantType.Insertion;
							break;
						case VariantType.SNV:
							variant.VarType = VariantType.SNVInsertion;
							break;
						case VariantType.Insertion:
							variant.VarType = VariantType.Insertion;
							break;
						case VariantType.Deletion:
							variant.VarType = VariantType.InsertionDeletion;
							break;
						default:
							variant.VarType = VariantType.Complex;
							break;
					}
					break;
				case VariantType.Deletion:
					switch (variant.VarType2)
					{
						case VariantType.Reference:
							variant.VarType = VariantType.Deletion;
							break;
						case VariantType.SNV:
							variant.VarType = VariantType.SNVDeletion;
							break;
						case VariantType.Insertion:
							variant.VarType = VariantType.InsertionDeletion;
							break;
						case VariantType.Deletion:
							variant.VarType = VariantType.Deletion;
							break;
						default:
							variant.VarType = VariantType.Complex;
							break;
					}
					break;
				default:
					variant.VarType = VariantType.Complex;
					break;
			}
		}

		/// <summary>
		/// Assign a variant type to a particular allele.  The rules are as follows:
		/// - If ref==alt, type is reference.  
		/// - Otherwise, trim off any common prefix and any common suffix.  Let |ref| denote the length of the
		///   reference allele after trimming, and |alt| denote the length of the alt allele after trimming.
		/// - If |ref|=0, it's an insertion
		/// - If |alt|=0, it's a deletion
		/// - If |ref|=|alt|=1, it's a SNV
		/// - If |ref| = |alt| > 1, it's a MNP
		/// - If |ref|>0 and |alt|>0 and |ref| != |alt|, it's a complex event
		/// </summary>
		private static VariantType GetAlleleVariantType(VcfVariant variant, int haplotype)
		{
			if (haplotype == 0)
				return VariantType.Reference;
			if (haplotype == -1)
				return VariantType.Missing;

			string altAllele = variant.VariantAlleles[haplotype - 1];
			return GetAlleleVariantType(variant.ReferenceAllele, altAllele);
		}

		/// <summary>
		/// Assign a variant type to a particular allele.  The rules are as follows:
		/// - If ref==alt, type is reference.  
		/// - Otherwise, trim off any common prefix and any common suffix.  Let |ref| denote the length of the
		///   reference allele after trimming, and |alt| denote the length of the alt allele after trimming.
		/// - If |ref|=0, it's an insertion
		/// - If |alt|=0, it's a deletion
		/// - If |ref|=|alt|=1, it's a SNV
		/// - If |ref| = |alt| > 1, it's a MNP
		/// - If |ref|>0 and |alt|>0 and |ref| != |alt|, it's a complex event
		/// </summary>
		public static VariantType GetAlleleVariantType(string referenceAllele, string altAllele)
		{
			if (referenceAllele == altAllele) return VariantType.Reference; // Sanity check; this should never trigger in practice
			if (altAllele == ".") return VariantType.Reference; // This shouldn't happen, but sometimes we see REF=A ALT=. GT=0/1
			if (referenceAllele.Length == 1 && altAllele.Length == 1) return VariantType.SNV;

			int commonPrefixBases = 0;
			for (; commonPrefixBases < referenceAllele.Length && commonPrefixBases < altAllele.Length; commonPrefixBases++)
			{
				if (referenceAllele[commonPrefixBases] != altAllele[commonPrefixBases]) break;
			}
			int commonSuffixBases = 0;
			for (; commonSuffixBases < referenceAllele.Length && commonSuffixBases < altAllele.Length; commonSuffixBases++)
			{
				if (referenceAllele[referenceAllele.Length - 1 - commonSuffixBases] !=
					altAllele[altAllele.Length - 1 - commonSuffixBases]) break;
			}
			int trimBases = commonPrefixBases + commonSuffixBases;
			if (trimBases > referenceAllele.Length) trimBases = referenceAllele.Length;
			if (trimBases > altAllele.Length) trimBases = altAllele.Length;
			if (referenceAllele.Length <= trimBases)
				return VariantType.Insertion;
			if (altAllele.Length <= trimBases)
				return VariantType.Deletion;
			if (trimBases + 1 == altAllele.Length && altAllele.Length == referenceAllele.Length)
				return VariantType.SNV;
			if (altAllele.Length == referenceAllele.Length)
				return VariantType.MNP;
			return VariantType.Complex;
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
			List<VcfVariant> allVariants = new List<VcfVariant>();
			using (VcfReader reader = new VcfReader(vcfPath))
			{
				foreach (VcfVariant variant in reader.GetVariants())
				{
					allVariants.Add(variant);
				}
			}
			return allVariants;
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