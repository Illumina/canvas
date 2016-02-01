
using System.Collections.Generic;

namespace Isas.Shared
{
    /// <summary>
    ///     Statistics specific to VeriSeqPGS
    /// </summary>
    // ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
    public class VeriSeqPGSStatistics : WorkflowStatistics
    {
        // ReSharper restore InconsistentNaming

        public override string FileName
        {
            get
            {
                return "VeriSeqRunStatistics.xml";
            }
        }

        public VeriSeqPGSStatistics()
        {
        }

        public VeriSeqPGSStatistics(int numberOfReads)
        {
            Version = 2;
            RunStats = new RunStatistics(numberOfReads);
        }
    }
}

