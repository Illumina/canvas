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
        public PhasedGenotype GenotypeTranspose()
        {
            return new PhasedGenotype(CopyNumberB, CopyNumberA);
        }

        /// <summary>
        /// Returns true when condition is satisfied and allele copy numbers are present
        /// </summary>
        /// <param name="phasedGenotype"></param>
        /// <returns></returns>
        public bool ContainsSharedAlleleA(PhasedGenotype phasedGenotype)
        {
            return CopyNumberA == phasedGenotype.CopyNumberA || CopyNumberA == phasedGenotype.CopyNumberB;
        }

        public bool ContainsSharedAlleleB(PhasedGenotype phasedGenotype)
        {
            return CopyNumberB == phasedGenotype.CopyNumberA || CopyNumberB == phasedGenotype.CopyNumberB;
        }

        public override int GetHashCode()
        {
            return CopyNumberA * 10 + CopyNumberB;
        }

        public override bool Equals(object obj)
        {
            return Equals(obj as PhasedGenotype);
        }

        private bool Equals(PhasedGenotype obj)
        {
            return obj != null && obj.CopyNumberA == CopyNumberA && obj.CopyNumberB == CopyNumberB; 
        }
    }
}