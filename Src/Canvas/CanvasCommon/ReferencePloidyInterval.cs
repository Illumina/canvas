using Isas.SequencingFiles;

namespace CanvasCommon
{
    public class ReferencePloidyInterval
    {
        public ReferenceInterval Interval { get; }
        public int ReferencePloidy { get; }

        public ReferencePloidyInterval(ReferenceInterval interval, int referencePloidy)
        {
            Interval = interval;
            ReferencePloidy = referencePloidy;
        }
    }
}