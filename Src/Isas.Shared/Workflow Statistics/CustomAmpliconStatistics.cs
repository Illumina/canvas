using System;
using System.Collections.Generic;
using System.Xml.Serialization;
using ProtoBuf;
namespace Isas.Shared
{
    /// <summary>
    ///     Statistics specific to the truseq amplicon workflow.
    /// </summary>
    // ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
    public class StatisticsAmplicon : WorkflowStatistics
    {
        public List<AmpliconStatistics> Amplicons; // across samples
        // this is big - will save this info in a separate, more compact file format
        [XmlIgnore]
        public List<List<AmpliconStatistics>> AmpliconsBySample;
        // ReSharper restore InconsistentNaming

        public override string WorkflowName
        {
            get { return "TruSeq Amplicon"; }
        }

        public StatisticsAmplicon()
        {
        }

        public StatisticsAmplicon(int numberOfReads)
        {
            RunStats = new RunStatistics(numberOfReads);
        }
    }

    /// <summary>
    ///     Used in Custom Amplicon workflow
    ///     We first generate one of these per sample + amplicon, then aggregate across samples into one per amplicon.
    /// </summary>
    // ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializable members that are sensitive to being changed
    [Serializable]
    [ProtoContract]
    public class AmpliconStatistics
    {
        [ProtoMember(1)]
        public long[] AlignedBaseCount; // array by read number
        [ProtoMember(2)]
        public long[] BaseCount; // array by read number
        [ProtoMember(3)]
        public long[] ClustersPF; // array by read number
        [ProtoMember(4)]
        public long[] ClustersRaw; // array by read number
        [ProtoMember(5)]
        public int ConcatenatedEnd;
        [ProtoMember(6)]
        public int ConcatenatedStart;
        [ProtoMember(7)]
        public string DesignChromosome;
        [ProtoMember(8)]
        public int DesignChromosomeStart;
        [ProtoMember(9)]
        public long[] ErrorCount; // array by read number
        [NonSerialized]
        public float[] ErrorRate; // array by read number
        [ProtoMember(10)]
        public string Manifest;
        [ProtoMember(11)]
        public string Name;
        [ProtoMember(12)]
        public long[] NoCallCount;
        [NonSerialized]
        public float[] NoCallRate; // array by read number
        [NonSerialized]
        public VariantsForSequence VariantInfo; // variant stats per amplicon
        // ReSharper restore InconsistentNaming

        // these coordinates are needed for the DataAccess interface
        // each amplicon has a start and end in the consolidated reference genome

        public AmpliconStatistics()
        {
        }

        public AmpliconStatistics(Dictionary<int, int> cyclesByRead)
        {
            int numberOfReads = cyclesByRead.Keys.Count;
            ClustersRaw = new long[numberOfReads];
            ClustersPF = new long[numberOfReads];
            ErrorRate = new float[numberOfReads];
            NoCallRate = new float[numberOfReads];
            ErrorCount = new long[numberOfReads];
            AlignedBaseCount = new long[numberOfReads];
            NoCallCount = new long[numberOfReads];
            BaseCount = new long[numberOfReads];
            VariantInfo = new VariantsForSequence();
        }

        public void AggregateFromProcess(AmpliconStatistics processAmplicon)
        {
            SampleStatistics.AggregateCounts(AlignedBaseCount, processAmplicon.AlignedBaseCount);
            SampleStatistics.AggregateCounts(BaseCount, processAmplicon.BaseCount);
            SampleStatistics.AggregateCounts(ClustersPF, processAmplicon.ClustersPF);
            SampleStatistics.AggregateCounts(ClustersRaw, processAmplicon.ClustersRaw);
            SampleStatistics.AggregateCounts(ErrorCount, processAmplicon.ErrorCount);
            SampleStatistics.AggregateCounts(NoCallCount, processAmplicon.NoCallCount);
        }

        /// <summary>
        /// </summary>
        /// <param name="ampliconName"></param>
        /// <returns></returns>
        public int GetTargetIndex(string ampliconName)
        {
            int idxLastDot = ampliconName.LastIndexOf('.');
            int ampliconIdx = 1;
            if (idxLastDot > 0)
                Int32.TryParse(ampliconName.Substring(idxLastDot + 1), out ampliconIdx);
            return ampliconIdx;
        }
    }
}