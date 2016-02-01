using System.Collections.Generic;
using System.Reflection;
using System.Xml.Serialization;

namespace Isas.Shared
{
	// ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
	public class RunStatistics
    {
        public string AnalysisSoftwareVersion;
        public DescriptiveStats[] ErrorRate; // by read, across all samples
        public DescriptiveStats[] NoCalls; // by read, across all samples
        public long NumberOfClustersPF;
        public long NumberOfClustersRaw;
        public long NumberOfDuplicateClusters;
        public long NumberOfUnalignedClusters;
        public long NumberOfUnalignedClustersPF;
        public long NumberOfUnindexedClusters;
        public long NumberOfUnindexedClustersPF;
        public float OverallCoverage; // across all samples and all reads
        public float[] OverallCoveragePerRead; // by read, across all samples

        public string OverallErrorPlotPath;
        public string OverallNoCallPlotPath;
        public ReadPairProperties OverallPairStatistics;

        [XmlIgnore] public Dictionary<string, List<EnrichmentStatistics.RegionStatistics>> RegionStatsBySample =
            new Dictionary<string, List<EnrichmentStatistics.RegionStatistics>>();

        public long YieldInBasesPF;
        public long YieldInBasesRaw;
		// ReSharper restore InconsistentNaming

        public RunStatistics()
        {
            string version = Assembly.GetExecutingAssembly().GetName().Version.ToString();
            AnalysisSoftwareVersion = version;
        }

        public RunStatistics(int numberOfReads)
        {
            ErrorRate = new DescriptiveStats[numberOfReads];
            NoCalls = new DescriptiveStats[numberOfReads];
            OverallCoveragePerRead = new float[numberOfReads];
            for (int readIndex = 0; readIndex < numberOfReads; readIndex++)
            {
                ErrorRate[readIndex] = new DescriptiveStats();
                NoCalls[readIndex] = new DescriptiveStats();
            }
            string version = Assembly.GetExecutingAssembly().GetName().Version.ToString();
            AnalysisSoftwareVersion = version;
        }
    }
}