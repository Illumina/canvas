namespace CanvasCommon
{
    public class Genotype
    {
        public int TotalCopyNumber { get; }

        public PhasedGenotype PhasedGenotype { get; }

        private Genotype(int totalCopyNumber, PhasedGenotype phasedGenotype)
        {
            TotalCopyNumber = totalCopyNumber;
            PhasedGenotype = phasedGenotype;
        }

        public static Genotype Create(int totalCopyNumber)
        {
            return new Genotype(totalCopyNumber, null);
        }

        public static Genotype Create(PhasedGenotype phasedGenotype)
        {
            return new Genotype(phasedGenotype.CopyNumberA + phasedGenotype.CopyNumberB, phasedGenotype);
        }

        public bool HasAlleleCopyNumbers => PhasedGenotype != null;

        public bool ContainsSharedAlleles(Genotype genotype)
        {
            return TotalCopyNumber == genotype.TotalCopyNumber;
        }

        public override int GetHashCode()
        {
            if (PhasedGenotype != null)
                return TotalCopyNumber * 100 + PhasedGenotype.GetHashCode();
            return TotalCopyNumber;
        }

        public override bool Equals(object obj)
        {
            return Equals(obj as Genotype);
        }
        private bool Equals(Genotype obj)
        {
            if (obj?.PhasedGenotype != null && PhasedGenotype != null)
                return obj.TotalCopyNumber == TotalCopyNumber && obj.PhasedGenotype.CopyNumberA == PhasedGenotype.CopyNumberA &&
                   obj.PhasedGenotype.CopyNumberB == PhasedGenotype.CopyNumberB;
            return obj != null && obj.TotalCopyNumber == TotalCopyNumber;
        }
    }
}