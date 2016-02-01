using System.Collections.Generic;
using System.Xml.Serialization;
using System;

namespace Isas.Shared
{
	// ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
	[Serializable]
    public class VariantStatistics
    {
        public string ChromosomeName = null; // if this is non-null then the statistics pertain to only this chromosome
        [XmlIgnore] public int EndPosition; // currently just used by pcr amplicon
        public string FileStub;
        public int NumberHeterozygousSNPs;
        public int NumberHomozygousSNPs;
        public int NumberOfDeletions;
        public int NumberOfIndels;
        public int NumberOfInsertions;
        public int NumberOfSNPs;
        public DescriptiveStats Qindel;
        public DescriptiveStats Qsnp;
        [XmlIgnore] public int StartPosition; // currently just used by pcr amplicon
        [XmlIgnore] public List<float> IndelQScores;
        [XmlIgnore] public List<float> SNVQScores;
		// ReSharper restore InconsistentNaming
	}

    public class SmallVariantStats
    {
        public VariantCountDetails SNVs;
        public VariantCountDetails Deletions;
        public VariantCountDetails Insertions;
        public ReadCountDetails ReadCounts;
	}

    //Variant Read Frequency is in the somatic snvs VCF file 
    //(it’s not the plot, that’s using the germline calls). 
    //It’s taking the tier1 counts of the ALT base divided by sum of tier1 counts of ALT + tier 1 counts of REF bases, in the Tumor sample column.
    public class ReadCountDetails
    {
        public int Tier1CountsOfAltSupport;
        public int Tier1CountsOfRefSupport;

        public int AltPlusRefSupport()
        {
            return (Tier1CountsOfAltSupport + Tier1CountsOfRefSupport);
        }

        //as percent
        public double AvgVariantAlleleFreqByReadCounts()
        {
            if (AltPlusRefSupport() == 0)
                return 0;

            return (100.0* (double)Tier1CountsOfAltSupport / (double)AltPlusRefSupport());
        }
    }

    public class VariantCountDetails
    {
        public int Total;
        public int Genic;
        public int Exonic;
        public int Coding;
        public int SpliceSite;
        public int StopGained;
        public int StopLost;
        public int Frameshift;
        public int NonSynonymous;
        public int Synonymous;
        public int MatureMicroRNA;
        public int UTR;
    }

    public class GenicCount
    {
		public int TotalVariants;
        public int TotalPassingVariants;
        public int GenicPassingVariants;
    }

    public class SVStats
    {
        public GenicCount Insertions;
        public GenicCount Deletions;
        public GenicCount Inversions;
        public GenicCount TandemDuplications;
        public GenicCount Breakends;
    }

    public class CNVStats
    {
        public GenicCount count;
    }

    /// <summary>
    ///     Utility class used in Amplicon workflow
    /// </summary>
	// ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
	public class VariantsForSequence
    {
        public List<string> Genotype = new List<string>();
        public List<string> Indel = new List<string>();
        public List<int> IndelPositions = new List<int>();
        public int NumberHeterozygousSNPs;
        public int NumberHomozygousSNPs;
        public int NumberOfDeletions;
        public int NumberOfIndels;
        public int NumberOfInsertions;
        public int NumberOfSNPs;
        public List<float> Qindel = new List<float>();
        public DescriptiveStats QindelStats;
        // for each SNP record position, GT and quality
        public List<float> Qsnp = new List<float>();
        public DescriptiveStats QsnpStats;
        public List<int> SNVPositions = new List<int>();
		// ReSharper restore InconsistentNaming

        // for each Indel record position, Indel and quality
    }
}