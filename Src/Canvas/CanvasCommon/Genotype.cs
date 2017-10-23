namespace CanvasCommon
{
    public class Genotype
    {
        public int CountsA { get; }
        public int CountsB { get; }
        
        public Genotype(int countsA, int countsB)
        {
            CountsA = countsA;
            CountsB = countsB;
        }
    }
}