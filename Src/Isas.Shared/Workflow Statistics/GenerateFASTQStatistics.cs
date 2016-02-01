using System.Collections.Generic;

namespace Isas.Shared
{
    /// <summary>
    ///     Statistics specific to 
    ///     STQ
    /// </summary>
	// ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
	public class StatisticsGenerateFASTQ : WorkflowStatistics
    {
		// ReSharper restore InconsistentNaming

		public StatisticsGenerateFASTQ()
        {
        }

        public StatisticsGenerateFASTQ(int numberOfReads)
        {
            Version = 2;
            RunStats = new RunStatistics(numberOfReads);
        }
    }
}
