using System.Collections.Generic;
using System.Xml.Serialization;
using System.Xml;
using ProtoBuf;
using InterOp.Model;
using InterOp.IO.Read;
using InterOp.IO.Read.Source;
using InterOp.IO.Utils;

// RunInfoRun and helper class RunInfoRead represent the run definition provided in the file RunInfo.xml from the run folder.
// Parsing of the .xml file is handled by the standard Illumina.InterOp libraries.
// The main information we need is the read structure (how many reads, how long, and which are index reads).  We also capture some
// other meta-data: Lane count, instrument name

namespace SequencingFiles
{
	// ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
	[ProtoContract]
	public class RunInfoRun
	{
		[ProtoMember(1)]
		public string RunID { get; set; }

		[ProtoMember(2)]
		public string Number { get; set; } // This attribute is available for P11, but not for historical MiSeq/HiSeq/GA runs

		[ProtoMember(3)]
		public string Flowcell { get; set; }

		[ProtoMember(4)]
		public string Instrument { get; set; }

		[ProtoMember(5)]
		public string Date { get; set; }

		[ProtoMember(6)]
		public List<RunInfoRead> Reads { get; set; }

		[ProtoMember(7)]
		public int LaneCount { get; set; }

		// ReSharper restore InconsistentNaming

		/// <summary>
		/// Deserialize RunInfo.xml using Illumina.InterOp, then parse the information into our objects.
		/// </summary>
		public static RunInfoRun Load(string runFolder)
		{
			var fileSystemSource = new FileSystemSource(runFolder);
			InterOp.IO.Read.Serialized.RunInfo run;
			using (var dataStream = fileSystemSource.GetRunInfo())
			{
				if (dataStream == null)
				{
					throw new System.Exception(string.Format("Error: RunInfo.xml file not found for run folder '{0}'", runFolder));
				}
				var serializer = new XmlSerializer(typeof(InterOp.IO.Read.Serialized.RunInfo));
				run = (InterOp.IO.Read.Serialized.RunInfo)serializer.Deserialize(dataStream);
			}
			RunInfoRun myRun = new RunInfoRun();
			myRun.RunID = run.Run.Id;
			if (run.Run.Number > 0)
				myRun.Number = run.Run.Number.ToString();

			myRun.Flowcell = run.Run.Flowcell;
			// Sanity check: Ensure flowcell isn't null!
			if (string.IsNullOrEmpty(myRun.Flowcell))
			{
				myRun.Flowcell = System.IO.Path.GetFileName(runFolder);
			}
			myRun.Instrument = run.Run.Instrument;
			myRun.Date = run.Run.Date;
			if (run.Run.FlowcellLayout != null)
				myRun.LaneCount = run.Run.FlowcellLayout.LaneCount;
			else
				myRun.LaneCount = 8;
			int totalCycles = 0;
			myRun.Reads = new List<RunInfoRead>();
			foreach (var read in run.Run.Reads)
			{
				RunInfoRead myRead = new RunInfoRead(read.NumCycles);
				myRead.IsIndex = read.IsIndexedRead == "Y" ? true : false;
				myRun.Reads.Add(myRead);
				totalCycles += read.NumCycles;
			}
			if (totalCycles <= 0)
				throw new System.Exception(string.Format("Error in RunInfo for {0} - no reads with 1 or more cycles found.  Please confirm this is a valid (version 2+) RunInfo.xml file", runFolder));
			return myRun;
		}

		/// <summary>
		/// remove me please! We only want smart RunInfos around here
		/// Once we completely switch to using SampleBcls as input and remove the AnalysisWorker.RunInfo dependency then this method can go away.
		/// </summary>
		public static RunInfoRun GetDummyRunInfo()
		{
			RunInfoRun run = new RunInfoRun();
			run = new RunInfoRun();
			run.Instrument = "UnknownInstrument";
			run.RunID = "UnknownID";
			run.Flowcell = "UnknownFlowcell";
			run.Reads = new List<RunInfoRead>();
			// Create dummy flow cell layout:
			run.LaneCount = 8;
			return run;
		}

		public int GetIndexReadNumber()
		{
			int readNumber = 1;
			foreach (RunInfoRead read in Reads)
			{
				if (read.IsIndex) return readNumber;
				readNumber++;
			}
			return 0;
		}

		public bool HasIndexRead()
		{
			foreach (RunInfoRead read in Reads)
			{
				if (read.IsIndex) return true;
			}
			return false;
		}

	}

	// ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
	[ProtoContract]
	public class RunInfoRead
	{
		[ProtoMember(1)]
		public int CycleCount { get; set; }

		[ProtoMember(2)]
		public int StartFromCycle { get; set; }

		[ProtoMember(3)]
		public int EndWithCycle { get; set; }

		public int UsedCycleCount { get { return System.Math.Max(0, EndWithCycle - StartFromCycle + 1); } }

		[ProtoMember(4)]
		public bool IsIndex { get; set; }
		// ReSharper restore InconsistentNaming

		public RunInfoRead(int cycles)
		{
			CycleCount = cycles;
			StartFromCycle = 1;
			EndWithCycle = cycles;
			IsIndex = false;
		}

		public RunInfoRead()
		{
		}

	}

}