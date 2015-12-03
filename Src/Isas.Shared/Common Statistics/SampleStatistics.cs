using System;
using System.Collections.Generic;
using System.Threading;
using ProtoBuf;

namespace Isas.Shared
{
	/// <summary>
	/// Serialized class for sample+chromosome stats.
	/// Resequencing workflow use one instance per SAMPLE + CHROMOSOME
	/// Amplicon workflow uses one instance per SAMPLE
	/// </summary>
	// ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
	[Serializable]
	[ProtoContract]
	public class SampleStatistics
	{
		[ProtoMember(1)]
		public float[] AverageAlignmentScore;
		[ProtoMember(2)]
		public float[] AverageErrorRate;
		[ProtoMember(3)]
		public float[] AverageNoCallRate;
		[ProtoMember(4)]
		public string Chromosome;
		[ProtoMember(5)]
		public long ClustersAlignedR1;
		[ProtoMember(6)]
		public long ClustersAlignedR2;
		[ProtoMember(7)]
		public string[] CoveragePlotPath;
		[ProtoMember(8)]
		public int[] E01;
		[ProtoMember(9)]
		public int[] E03;
		[ProtoMember(10)]
		public int[] E05;
		[ProtoMember(11)]
		public int[] E10;
		[ProtoMember(12)]
		public string[] ErrorCoveragePlotPath;
		[ProtoMember(13)]
		public string ErrorPlotPath;
		[ProtoMember(14)]
		public float[] GCBias;
		[ProtoMember(15)]
		public string[] GCBiasPlotPath;
		[ProtoMember(16)]
		public string GenomeName;
		[ProtoMember(17)]
		public int GenomeNumber;
		[ProtoMember(18)]
		public string NoCallPlotPath;
		[ProtoMember(19)]
		public long NumberOfClustersPF;
		[ProtoMember(20)]
		public long NumberOfClustersRaw;
		[ProtoMember(21)]
		public long[] PerfectReads;
		[ProtoMember(22)]
		public DescriptiveStats SampleCoverage;
		[ProtoMember(23)]
		public DescriptiveStats[] SampleCoverageByRead;
		[ProtoMember(24)]
		public string SampleID; // user-specified
		[ProtoMember(25)]
		public string SampleName; // user-specified
		[ProtoMember(26)]
		public int SampleNumber;
		[NonSerialized]
		public VariantStatistics VariantStats;
		// ReSharper restore InconsistentNaming

		public SampleStatistics()
		{
		}

		public SampleStatistics(int numberOfReads)
		{
			if (numberOfReads > 0)
			{
				E01 = new int[numberOfReads];
				E03 = new int[numberOfReads];
				E05 = new int[numberOfReads];
				E10 = new int[numberOfReads];
				PerfectReads = new long[numberOfReads];
				AverageErrorRate = new float[numberOfReads];
				AverageNoCallRate = new float[numberOfReads];
				GCBias = new float[numberOfReads];
				GCBiasPlotPath = new string[numberOfReads];
				CoveragePlotPath = new string[numberOfReads];
				ErrorCoveragePlotPath = new string[numberOfReads];
				SampleCoverageByRead = new DescriptiveStats[numberOfReads];
				for (int index = 0; index < numberOfReads; index++)
				{
					SampleCoverageByRead[index] = new DescriptiveStats();
				}
				AverageAlignmentScore = new float[numberOfReads];
			}
		}

		/// <summary>
		///     Helper function to aggregate statistics when counting was split among multiple threads/processes
		/// </summary>
		public void AggregateFromProcess(SampleStatistics sampleStats)
		{
			Interlocked.Add(ref ClustersAlignedR1, sampleStats.ClustersAlignedR1);
			Interlocked.Add(ref ClustersAlignedR2, sampleStats.ClustersAlignedR2);
			AggregateCounts(PerfectReads, sampleStats.PerfectReads);
		}

		/// <summary>
		///     Helper function to aggregate counts when counting was split among multiple threads/processes
		/// </summary>
		public static void AggregateCounts(float[] aggregate, float[] part)
		{
			if (aggregate == null || part == null)
				return;

			// cannot do Interlocked.Increment on floats so just lock down the whole array
			// this should be the only function that modifies the array from multiple threads
			lock (aggregate)
			{
				for (int i = 0; i < aggregate.Length; i++)
				{
					aggregate[i] += part[i];
				}
			}
		}

		/// <summary>
		///     Helper function to aggregate counts when counting was split among multiple threads/processes
		/// </summary>
		public static void AggregateCounts(long[] aggregate, long[] part)
		{
			if (aggregate == null || part == null)
				return;

			for (int i = 0; i < aggregate.Length; i++)
			{
				Interlocked.Add(ref aggregate[i], part[i]);
			}
		}

		/// <summary>
		///     Helper function to aggregate counts when counting was split among multiple threads/processes
		/// </summary>
		public static void AggregateCounts(IDictionary<long, long> aggregate, IDictionary<long, long> part)
		{
			if (aggregate == null || part == null)
				return;

			lock (aggregate)
			{
				foreach (KeyValuePair<long, long> kvp in part)
				{
					long value;
					if (aggregate.TryGetValue(kvp.Key, out value))
						aggregate[kvp.Key] = value + kvp.Value;
					else
						aggregate[kvp.Key] = kvp.Value;
				}
			}
		}

		/// <summary>
		///     Helper function to aggregate counts when counting was split among multiple threads/processes
		/// </summary>
		public static void AggregateCounts(int[] aggregate, int[] part)
		{
			if (aggregate == null || part == null)
				return;

			for (int i = 0; i < aggregate.Length; i++)
			{
				Interlocked.Add(ref aggregate[i], part[i]);
			}
		}
	}

	/// <summary>
	/// Serialized class - overall stats per sample (rather than per sample+chromosome)
	/// </summary>
	// ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
	public class SummarizedSampleStatistics
	{
		public long ClustersAlignedR1;
		public long ClustersAlignedR2;
		public float ErrorRateR1;
		public float ErrorRateR2;
		public int FragmentLengthMax;
		public int FragmentLengthMedian;
		public int FragmentLengthMin;
		public uint FragmentLengthSD;
		public float NoCallRateR1;
		public float NoCallRateR2;
		public int NumberHeterozygousSNPs;
		public int NumberHomozygousSNPs;
		public long NumberOfClustersPF;
		public long NumberOfClustersRaw;
		public int NumberOfDeletions;
		public int NumberOfInsertions;
		public string SampleID;
		public string SampleName;
		public int SampleNumber;
		public VariantConcordance Variants;
		public float WeightedCoverage;
		public float PercentQ30; // %Q30 for all PF reads
		public float PercentQ30R1;
		public float PercentQ30R2;
		public float PercentQ30Aligned; // %Q30 for all aligned reads
		public float PercentQ30AlignedR1;
		public float PercentQ30AlignedR2;
		public long[] QHistogramR1;
		public long[] QHistogramR2;
		public long[] QHistogramAlignedR1;
		public long[] QHistogramAlignedR2;
		public SmallVariantStats SmallVariants;
		public SVStats StructuralVariants;
		public CNVStats CopyNumberVariants;

		// PUMA metrics
        // Percent callability is computed for non decoy chromosomes
		public float? PercentCallability;
		public float? PercentContamination;
		public float? PercentArrayConcordance;

		// autosome
		public float? MeanAutosomeCoverage;
		public float? PercentAutosomeBasesAtLeast1x;
		public float? PercentAutosomeBasesAtLeast10x;
		public float? PercentAutosomeBasesAtLeast15x;
		public float? PercentAutosomeCallability;

		// autosome exon
		public float? MeanAutosomeExonCoverage;
		public float? PercentAutosomeExonBasesAtLeast1x;
		public float? PercentAutosomeExonBasesAtLeast10x;
		public float? PercentAutosomeExonBasesAtLeast15x;
		public float? PercentAutosomeExonCallability;
		// ReSharper restore InconsistentNaming

		//add these methods so we don't serialize values that are not provided
		public bool ShouldSerializePercentCallability()
		{
			return PercentCallability.HasValue;
		}
		public bool ShouldSerializePercentContamination()
		{
			return PercentContamination.HasValue;
		}
		public bool ShouldSerializePercentArrayConcordance()
		{
			return PercentArrayConcordance.HasValue;
		}
		public bool ShouldSerializeMeanAutosomeCoverage()
		{
			return MeanAutosomeCoverage.HasValue;
		}
		public bool ShouldSerializePercentAutosomeBasesAtLeast1x()
		{
			return PercentAutosomeBasesAtLeast1x.HasValue;
		}
		public bool ShouldSerializePercentAutosomeBasesAtLeast10x()
		{
			return PercentAutosomeBasesAtLeast10x.HasValue;
		}
		public bool ShouldSerializePercentAutosomeBasesAtLeast15x()
		{
			return PercentAutosomeBasesAtLeast15x.HasValue;
		}
		public bool ShouldSerializeMeanAutosomeExonCoverage()
		{
			return MeanAutosomeExonCoverage.HasValue;
		}
		public bool ShouldSerializePercentAutosomeExonBasesAtLeast1x()
		{
			return PercentAutosomeExonBasesAtLeast1x.HasValue;
		}
		public bool ShouldSerializePercentAutosomeExonBasesAtLeast10x()
		{
			return PercentAutosomeExonBasesAtLeast10x.HasValue;
		}
		public bool ShouldSerializePercentAutosomeExonBasesAtLeast15x()
		{
			return PercentAutosomeExonBasesAtLeast15x.HasValue;
		}
		public bool ShouldSerializePercentAutosomeCallability()
		{
			return PercentAutosomeCallability.HasValue;
		}
		public bool ShouldSerializePercentAutosomeExonCallability()
		{
			return PercentAutosomeExonCallability.HasValue;
		}
	}

	/// <summary>
	///     keeps dbSNP concordance and ti/tv ratio until statistics evaluation
	/// </summary>
	// ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
	public class VariantConcordance
	{
		public VariantTypeStatistics Deletions;
		public VariantTypeStatistics Insertions;

		public int TransitionCount;
		public int TransversionCount;
		public VariantTypeStatistics SNPs;
		public double TransitionTransversionRatio;
		// ReSharper restore InconsistentNaming

		// constructor
		public VariantConcordance()
		{
			SNPs = new VariantTypeStatistics();
			Insertions = new VariantTypeStatistics();
			Deletions = new VariantTypeStatistics();
		}

		/// <summary>
		///     calculates the derived statistics (ti/tv ratio, % dbSNP concordance, etc.)
		/// </summary>
		public void Calculate()
		{
			SNPs.Calculate();
			Insertions.Calculate();
			Deletions.Calculate();
			TransitionTransversionRatio = (TransversionCount > 0 ? TransitionCount / (double)TransversionCount : 0.0);
		}

		// ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
		public class VariantTypeStatistics
		{
			public double HeterozygoteToHomozygoteRatio;
			public int NumHeterozygousVariants;
			public int NumHomozygousVariants;
			public int NumTotalVariants; //passing and non-passing
			public int NumPassingVariants;
			public int NumVariantsInDbSnp;
			public int NumVariantsNotInDbSnp;
			public double PercentVariantsInDbSnp;
			public double PercentVariantsNotInDbSnp;
			// ReSharper restore InconsistentNaming

			// constructor

			/// <summary>
			///     calculates the % dbSNP concordance
			/// </summary>
			public void Calculate()
			{
				double totalChecked = NumVariantsInDbSnp + NumVariantsNotInDbSnp;
				PercentVariantsInDbSnp = (totalChecked > 0 ? NumVariantsInDbSnp / totalChecked * 100.0 : 0.0);
				PercentVariantsNotInDbSnp = (totalChecked > 0 ? NumVariantsNotInDbSnp / totalChecked * 100.0 : 0.0);
				HeterozygoteToHomozygoteRatio = (NumHomozygousVariants > 0
													 ? NumHeterozygousVariants / (double)NumHomozygousVariants
													 : 0.0);
			}
		}
	}
}