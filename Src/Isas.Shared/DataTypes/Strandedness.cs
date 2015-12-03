namespace Isas.Shared.DataTypes
{
    /// <summary>
    /// Represents the orientation of the RNA in relation to the transcript strand
    /// </summary>
    public enum Strandedness
    {
        Unstranded, First, Second
    }

    public static class StrandedExtensions
    {
        public static bool IsStranded(this Strandedness stranded)
        {
            return stranded != Strandedness.Unstranded;
        }
    }
}
