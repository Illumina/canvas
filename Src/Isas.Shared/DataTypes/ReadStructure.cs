using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using System.Threading.Tasks;

namespace Illumina.SecondaryAnalysis
{
	public class ReadStructure
	{
		public readonly List<Read> Reads;
		public readonly int? Read1StartFromCycle;
		public readonly int? Read1EndWithCycle;
		public readonly int? Read2StartFromCycle;
		public readonly int? Read2EndWithCycle;

		public Read Read1
		{
			get
			{
				return Reads.FirstOrDefault(read => read == GetFirstMax() || read == GetSecondMax());
			}
		}

		public Read Index1
		{
			get
			{
				return Reads.FirstOrDefault(r => r.IsIndex);
			}
		}

		public Read Read2
		{
			get
			{
				return Reads.Where(read => read == GetFirstMax() || read == GetSecondMax())
							.Skip(1)
							.FirstOrDefault();
			}
		}

		public Read Index2
		{
			get
			{
				return Reads.Where(r => r.IsIndex)
							.Skip(1)
							.FirstOrDefault();
			}
		}

		private Read GetFirstMax()
		{
			return Reads.Where(read => !read.IsIndex)
						.OrderByDescending(read => read.Length)
						.FirstOrDefault();
		}

		private Read GetSecondMax()
		{
			return Reads.Where(read => !read.IsIndex)
						.OrderByDescending(read => read.Length)
						.Skip(1)
						.FirstOrDefault();
		}

		public ReadStructure(List<Read> reads, int? read1StartFromCycle, int? read1EndWithCycle, int? read2StartFromCycle, int? read2EndWithCycle)
		{
			Reads = reads;
			Read1StartFromCycle = read1StartFromCycle;
			Read1EndWithCycle = read1EndWithCycle;
			Read2StartFromCycle = read2StartFromCycle;
			Read2EndWithCycle = read2EndWithCycle;
			if (Reads.Any(r => r == null))
				throw new ArgumentException("Reads cannot be null");
			if (Reads.Count(r => r.IsIndex) > 2)
				throw new ArgumentException("A maximum of 2 indexed reads is supported");
			if (Read1 == null)
				throw new ArgumentException("At least 1 non-indexed read is required");
		}

		public class Read
		{
			public readonly int Length;
			public readonly bool IsIndex;

			public Read(int length, bool isIndex)
			{
				Length = length;
				IsIndex = isIndex;
			}
		}
	}
}
