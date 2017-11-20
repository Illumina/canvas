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
        /// <summary>
        /// Returns true when condition is satisfied and allele copy numbers are present
        /// </summary>
        /// <param name="phasedGenotype"></param>
        /// <returns></returns>
        public bool ContainsSharedAlleles(PhasedGenotype phasedGenotype)
        {
            return CopyNumberA == phasedGenotype.CopyNumberA | CopyNumberA == phasedGenotype.CopyNumberB | CopyNumberB == phasedGenotype.CopyNumberA | CopyNumberB == phasedGenotype.CopyNumberB;
        }

    }
}