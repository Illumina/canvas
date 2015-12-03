using System;
using System.Collections.Generic;
using System.Xml.Serialization;
using ProtoBuf;

namespace Isas.Shared
{
    /// <summary>
    ///     Statistics specific to PCR Amplicon
    /// </summary>
    // ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
    public class PCRAmpliconStatistics : WorkflowStatistics
    {
        public List<RegionStatistics> Regions = new List<RegionStatistics>();
        // ReSharper restore InconsistentNaming

        public override string FileName
        {
            get { return "PCRAmpliconStatistics.xml"; }
        }

        public PCRAmpliconStatistics()
        {
        }

        public PCRAmpliconStatistics(int numberOfReads)
        {
            RunStats = new RunStatistics(numberOfReads);
        }
    }

    /// <summary>
    ///     Used in Amplicon workflow
    ///     One of these per sample
    /// </summary>
    // ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
    [Serializable]
    [ProtoContract]
    public class RegionStatistics
    {
        [XmlIgnore]
        [ProtoMember(1)]
        public long[] AlignedBaseCount;
        [XmlIgnore]
        [ProtoMember(2)]
        public int[] BaseCount;
        [ProtoMember(3)]
        public string Chromosome;
        [ProtoMember(4)]
        public long[] ClustersPF;
        [ProtoMember(5)]
        public long[] ClustersRaw;
        [ProtoMember(6)]
        public float Coverage;
        [XmlIgnore]
        [ProtoMember(7)]
        public int[] CoverageByPosition;
        [ProtoMember(8)]
        public int EndPosition;
        [XmlIgnore]
        [ProtoMember(9)]
        public long[] ErrorCount;
        [ProtoMember(10)]
        public float[] ErrorRate;
        [ProtoMember(11)]
        public int Genome;
        [ProtoMember(12)]
        public string Manifest;
        [XmlIgnore]
        [ProtoMember(13)]
        public int[] NoCallCount;
        [ProtoMember(14)]
        public float[] NoCallRate;

        [XmlIgnore]
        [ProtoMember(15)]
        public Dictionary<string, int> ReadsBySample; // Aligned reads by sample - this counts R1 and R2 separately, so this is up to 2x the sample cluster count
        [ProtoMember(16)]
        public string RegionName;

        // these coordinates are needed for the DataAccess interface
        // each amplicon has a start and end in the consolidated reference genome
        [ProtoMember(17)]
        public int StartPosition;
        [NonSerialized]
        public VariantsForSequence VariantInfo;
        // ReSharper restore InconsistentNaming

        public RegionStatistics()
        {
        }

        //for cloning
        public RegionStatistics(RegionStatistics regionStats)
        {
            AlignedBaseCount = regionStats.AlignedBaseCount;
            BaseCount = regionStats.BaseCount;
            Chromosome = regionStats.Chromosome;
            ClustersPF = regionStats.ClustersPF;
            ClustersRaw = regionStats.ClustersRaw;
            Coverage = regionStats.Coverage;
            CoverageByPosition = regionStats.CoverageByPosition;
            EndPosition = regionStats.EndPosition;
            ErrorCount = regionStats.ErrorCount;
            ErrorRate = regionStats.ErrorRate;
            Genome = regionStats.Genome;
            Manifest = regionStats.Manifest;
            NoCallCount = regionStats.NoCallCount;
            NoCallRate = regionStats.NoCallRate;
            ReadsBySample = regionStats.ReadsBySample;
            RegionName = regionStats.RegionName;
            StartPosition = regionStats.StartPosition;
            VariantInfo = regionStats.VariantInfo;
        }

        public RegionStatistics(Dictionary<int, int> cyclesByRead, string chromosome, string manifest, int genome,
                                string name, int start, int end)
        {
            Chromosome = chromosome;
            Manifest = manifest;
            StartPosition = start;
            EndPosition = end;
            RegionName = name;
            Genome = genome;

            int numberOfReads = cyclesByRead.Keys.Count;
            ClustersRaw = new long[numberOfReads];
            ClustersPF = new long[numberOfReads];
            ErrorRate = new float[numberOfReads];
            NoCallRate = new float[numberOfReads];
            ErrorCount = new long[numberOfReads];
            NoCallCount = new int[numberOfReads];
            AlignedBaseCount = new long[numberOfReads];
            BaseCount = new int[numberOfReads];
            CoverageByPosition = new int[EndPosition - StartPosition + 1];
            ReadsBySample = new Dictionary<string, int>();
            VariantInfo = new VariantsForSequence();
        }

        public string GetKey()
        {
            return string.Format("{0}.{1}", Manifest, RegionName);
        }
    }
}