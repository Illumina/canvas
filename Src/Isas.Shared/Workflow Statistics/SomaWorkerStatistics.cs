using System.Collections.Generic;
using System.IO;
using System.Xml.Serialization;

namespace Isas.Shared
{
    /// <summary>
    ///     Statistics specific to SomaWorker
    /// </summary>
    // ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed   
    public class StatisticsSomaWorker : WorkflowStatistics
    {
        // ReSharper restore InconsistentNaming

        public override string FileName
        {
            get { return "AmpliconRunStatistics.xml"; }
        }

        public override string WorkflowName
        {
            get { return "Amplicon - DS"; }
        }

        public StatisticsSomaWorker()
        {
        }

        public StatisticsSomaWorker(int numberOfReads)
        {
            RunStats = new RunStatistics(numberOfReads);
        }
    }
}

