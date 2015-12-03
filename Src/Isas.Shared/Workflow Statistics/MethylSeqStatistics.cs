
namespace Isas.Shared
{
	// ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
	public class MethylSeqStatistics : WorkflowStatistics
    {
		// ReSharper restore InconsistentNaming

		public MethylSeqStatistics()
        {
        }

        public MethylSeqStatistics(int numberOfReads)
        {
            RunStats = new RunStatistics(numberOfReads);
        }
    }
}
