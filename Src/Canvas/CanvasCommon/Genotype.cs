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
    }
}