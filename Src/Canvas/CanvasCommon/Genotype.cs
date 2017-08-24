namespace CanvasCommon
{
    public class Genotype
    {
        public int CountsA { get; }
        public int CountsB { get; }

        public Genotype()
        {
            CountsA = 0;
            CountsB = 0;
        }

        public Genotype(int countsA, int countsB)
        {
            CountsA = countsA;
            CountsB = countsB;
        }
    }
}