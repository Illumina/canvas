using System;
using System.IO;
using System.Linq;
using System.Xml;
using System.Xml.Serialization;
using SequencingFiles;
using System.Collections.Generic;


namespace Isas.Shared
{
    /// <summary>
    /// Simple wrapper for WorkflowStatisticsCore which adds in RunStatistics.  
    /// CompletedJobInfo.xml is a serialized AnalysisJobInfo which inherits from WorkflowStatisticsCore
    /// ResequencingRunStatistics.xml (and other similar files) are serialized versions of files which inherit from WorkflowStatistics.
    /// </summary>
    public abstract class WorkflowStatistics : WorkflowStatisticsCore
    {
        public RunStatistics RunStats = new RunStatistics();

        public List<SummarizedSampleStatistics> OverallSamples = new List<SummarizedSampleStatistics>(); // summarizedSampleStatisticsList
        public List<ReadPairProperties> PairedEndByGenome = new List<ReadPairProperties>(); // by genome
        public List<SampleStatistics> StatsSamples = new List<SampleStatistics>();    //it is sampleChromosomeStatisticsList

        public virtual string FileName
        {
            get { return GetType().Name.Replace("Statistics", "") + "RunStatistics.xml"; }
        }

        public virtual string WorkflowName
        {
            get { return GetType().Name.Replace("Statistics", ""); }
        }

        /// <summary>
        /// Get list of all the SampleStatistics for a given sample number
        /// Note this is not very efficient way to access this information. To make it faster implement a dictionary lookup inside this class
        /// </summary>
        public List<SampleStatistics> GetSampleChromosomeStatisticsListOfOneSample(int sampleNumber)
        {
            List<SampleStatistics> sampleChromosomeStatisticsList = new List<SampleStatistics>();
            foreach (SampleStatistics sampleChromosomeStatistics in this.StatsSamples)
            {
                if (sampleChromosomeStatistics.SampleNumber == sampleNumber)
                    sampleChromosomeStatisticsList.Add(sampleChromosomeStatistics);
            }
            return sampleChromosomeStatisticsList;
        }

        public void CleanDummyChromosomeStatsRecords()
        {
            StatsSamples.RemoveAll(x => string.IsNullOrEmpty(x.Chromosome));
        }

        public void CopySampleChromosomeStatisticsToSummarizedSampleStatistics()
        {

        }

    }

    // ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
    public abstract class WorkflowStatisticsCore
    {
        #region Serializeable Types and Members
        public string AnalysisFolder; // the alignment folder for this job
        public DateTime CompletionTime;
        public string Error { get; set; } // if this is non-null then the job was terminated for some reason and didn't complete
        public string Warning { get; set; }
        public RunInfoRun RTARunInfo; // TODO: This is redundant with AnalysisWorker.RTARunInfo. What should go here if more than one runFolder is used. Todo: Eventually stop using RNARunInfo completely
        public string RunFolder; // the run folder
		[XmlIgnore]
		public List<SampleSheet.Sample> Samples; // TODO: make this the only copy in memory
	    public string WorkflowType;
        public DateTime StartTime;

        [XmlAttribute]
        public int Version;
        public SecondaryAnalysisWorkflow Workflow;

        #endregion
        // ReSharper restore InconsistentNaming


    }

}