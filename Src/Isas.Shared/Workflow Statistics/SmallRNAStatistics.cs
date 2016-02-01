using System.Collections.Generic;

namespace Isas.Shared
{
    // ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
    public class SmallRNAStatistics
    {
        public long ClusterCount;
        public long ClustersAlignedAbundant;
        public long ClustersAlignedGenome;
        public long ClustersAlignedMiRNA;
        public long ClustersAlignedRNA;
        public long ClustersPF;
        public long ClustersUnaligned;
        public string Sample;
        public string SampleID;
        public int SampleNumber;
        // ReSharper restore InconsistentNaming
    }

    public class SmallRNAReadCountStats
    {
        public int ClusterCount;
        public int ClusterPFCount;
    }

    // ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
    public class StatisticsSmallRNA : WorkflowStatistics
    {
        public List<SmallRNAStatistics> StatisticsBySample;
        // ReSharper restore InconsistentNaming

        public StatisticsSmallRNA()
        {
            //these are not used so set them to null so that they aren't serialized
            OverallSamples = null;
            PairedEndByGenome = null; // by genome
            StatsSamples = null;
        }
    }
}