using System.Collections.Generic;

namespace Isas.Shared
{
    /// <summary>
    ///     Statistics specific to RNAQuant
    /// </summary>
    // ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
    public class StatisticsRNAQuant : WorkflowStatistics
    {
        // ReSharper restore InconsistentNaming

        public override string WorkflowName
        {
            get
            {
                return "RNAQuantification";
            }
        }
    }
}
