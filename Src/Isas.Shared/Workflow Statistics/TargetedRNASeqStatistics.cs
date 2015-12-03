using System.Collections.Generic;
using System.Xml.Serialization;

namespace Isas.Shared
{
    /// <summary>
    ///     Statistics specific to the custom amplicon workflow.
    /// </summary>
    // ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
    public class StatisticsTargetedRNASeq : WorkflowStatistics
    {
        public List<AmpliconStatistics> Amplicons; // across samples
        // this is big - will save this info in a separate, more compact file format
        [XmlIgnore]
        public List<List<AmpliconStatistics>> AmpliconsBySample;
        // ReSharper restore InconsistentNaming

        public override string WorkflowName
        {
            get
            {
                return "TargetedRnaSeq";
            }
        }

        public StatisticsTargetedRNASeq()
        {
        }

        public StatisticsTargetedRNASeq(int numberOfReads)
        {
            RunStats = new RunStatistics(numberOfReads);
        }
    }
}