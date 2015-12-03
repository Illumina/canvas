using Isas.Shared;
using Newtonsoft.Json;
using SequencingFiles;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Illumina.SecondaryAnalysis
{
    public class Fastq : IMoveable<Fastq>, IMoveableResult<Fastq>
	{
		public readonly IFileLocation Read1;
		public readonly IFileLocation Read2;

        public long? ReadCount
        {
            get; private set;
        }
		public int? Read1Length
        {
            get; private set;
        }
        public int? Read2Length
        {
            get; private set;
        }

        public readonly Adapters Adapters; // In case the Fastq haven't been trimmed yet
        public readonly bool IsSanitized;

        [JsonIgnore]
        public bool IsPairedEnd => Read2 != null;
		
		public Fastq(IFileLocation read1, IFileLocation read2 = null, long? readCount = null, int? read1Length = null,
			int? read2Length = null, Adapters adapters = null, bool sanitized = true)
		{
			if (read1 == null)
				throw new ArgumentException("Fastq read 1 cannot be null");
			if (read2 == null && read2Length.HasValue)
				throw new ArgumentException("Read length 2 cannot be specified when there is no read 2 fastq file");

			Read1 = read1;
			Read2 = read2;
            Read1Length = read1Length;
            Read2Length = read2Length;
            
			ReadCount = readCount;

			Adapters = adapters;
            IsSanitized = sanitized;
		}

		public long GetReadCount()
		{
			int readCount = 0;
			using (FastqReader reader = new FastqReader(Read1.FullName))
			{
				BoltRead read = new BoltRead();
				while (reader.GetNextFastqEntry(ref read)) readCount++;
			}
			return readCount;
		}

        public int ComputeRead1Length()
        {
            return GetReadLength(Read1).Value;
        }
        public int? ComputeRead2Length()
        {
            return GetReadLength(Read2);
        }

        public static int? GetReadLength(IFileLocation path, int num_reads=10000)
		{
			if (path == null) return null;
			int readLength = 0;
			using (FastqReader reader = new FastqReader(path.FullName))
			{
				BoltRead read = new BoltRead();
				for (int readNumber = 0; readNumber < num_reads; readNumber++)
				{
					if (!reader.GetNextFastqEntry(ref read)) break;
					readLength = Math.Max(readLength, read.Bases.Length);
				}
			}
			return readLength;
		}

		/// <summary>
		/// Equal to Read1Count for paired-end data, null for single-end data
		/// </summary>
		[JsonIgnore]
		public long? Read2Count
		{
			get
			{
				if (Read2 == null) return null;
				return ReadCount;
			}
		}
        
		/// <summary>
		/// Upper limit on the number of bases contained in these fastq files. Includes no-calls (N). 
		/// Includes adapter-trimmed bases which were removed from the fastq files.
		/// </summary>
		[JsonIgnore]
		public long? MaxBaseCount
		{
			get
			{
                long? maxBaseCount = null;
                if (ReadCount.HasValue && Read1Length.HasValue)
				    maxBaseCount = ReadCount.Value * Read1Length.Value;
                if (IsPairedEnd)
                    if (Read2Length.HasValue)
                        maxBaseCount += ReadCount.Value * Read2Length.Value;
                    else
                        maxBaseCount *= 2;
				return maxBaseCount;
			}
		}

		public Fastq Move(FileNamingConvention getDestination)
		{
			IFileLocation newRead1 = Read1.Move(getDestination);

			IFileLocation newRead2;
			if (Read2 != null)
			{
				newRead2 = Read2.Move(getDestination);
			}
			else { newRead2 = null; }

			return new Fastq(newRead1, newRead2, ReadCount, Read1Length, Read2Length, Adapters);
		}

		public Fastq MoveOld(IDirectoryLocation directory)
		{
			return Move(directory.GetFileLocation);
		}

		public delegate IFileLocation GetOutputPath(string currentFileName);
		public Fastq MoveOld(GetOutputPath getOutputPath)
		{
			var outputRead1 = getOutputPath(Read1.Name);
			Read1.MoveTo(outputRead1);
			IFileLocation outputRead2 = null;
			if (Read2 != null)
			{
				outputRead2 = getOutputPath(Read2.Name);
				Read2.MoveTo(outputRead2);
			}
			return new Fastq(outputRead1, outputRead2, ReadCount, Read1Length, Read2Length, Adapters);
		}

        public Fastq Move(Fastq destination)
        {
            // Move files, return new Fastq
            IFileLocation newRead1 = Read1.MoveAndLink(destination.Read1);

            IFileLocation newRead2;
            if (Read2 != null)
            {
                newRead2 = Read2.MoveAndLink(destination.Read2);
            }
            else { newRead2 = null; }

            return new Fastq(newRead1, newRead2, ReadCount, Read1Length, Read2Length, Adapters);
        }
    }

	public static class FastqMover
	{
		public delegate IFileLocation GetOutputPath(SampleInfo info, int fastqIndex, string currentFileName);

		/// <summary>
		/// Move sample-specific collections of fastq files to a new location as defined by a delegate naming convention
		/// </summary>
		/// <param name="sampleSet">The sample-specific collection of fastq files to move</param>
		/// <param name="getOutputPath">A delegate for getting the final path for individual fastq files after the move</param>
		/// <returns>The sample-specific collection of fastq files in their new location after the move</returns>
		public static SampleSet<IEnumerable<Fastq>> MoveOld(this SampleSet<IEnumerable<Fastq>> sampleSet, GetOutputPath getOutputPath)
		{
			//use ToList() here to prevent lazy evaluation which would delay the Move!
			return sampleSet.SelectData(
				(info, fastqs) =>
					(IEnumerable<Fastq>)fastqs.Select(
						(fastq, index) =>
							fastq.MoveOld(currentFileName => getOutputPath(info, index, currentFileName))).ToList());
		}
	}
}
