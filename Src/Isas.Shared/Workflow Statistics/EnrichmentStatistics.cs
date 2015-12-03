namespace Isas.Shared
{
    /// <summary>
    ///     Statistics specific to Enrichment workflow
    /// </summary>
    // ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
    public class EnrichmentStatistics : PCRAmpliconStatistics
    {
        // ReSharper restore InconsistentNaming
        public override string FileName { get { return "EnrichmentStatistics.xml"; } }

        // constructor
        public EnrichmentStatistics()
        {
        }

        // constructor
        public EnrichmentStatistics(int numberOfReads)
            : this()
        {
            RunStats = new RunStatistics(numberOfReads);
        }
    }
}