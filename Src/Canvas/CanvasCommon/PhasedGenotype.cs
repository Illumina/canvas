using Illumina.Common;

namespace CanvasCommon
{
    public class PhasedGenotype
    {
        public int CopyNumberA { get; }
        public int CopyNumberB { get; }

        public PhasedGenotype(int copyNumberA, int copyNumberB)
        {
            CopyNumberA = copyNumberA;
            CopyNumberB = copyNumberB;
        }
    }
}