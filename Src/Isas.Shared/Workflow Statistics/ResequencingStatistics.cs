using System.Collections.Generic;

namespace Isas.Shared
{
    /// <summary>
    ///     Statistics specific to Resequencing
    /// </summary>
    // ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
    public class StatisticsResequencing : WorkflowStatistics
    {
        // ReSharper restore InconsistentNaming

        public StatisticsResequencing()
        {
        }

        public StatisticsResequencing(int numberOfReads)
        {
            Version = 2;
            RunStats = new RunStatistics(numberOfReads);
        }
    }
}