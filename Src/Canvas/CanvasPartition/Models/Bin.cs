namespace CanvasPartition.Models
{
    public class Bin
    {
        public uint Start;
        public uint End;
        public double Coverage;

        public Bin(uint start, uint end, double coverage)
        {
            Start = start;
            End = end;
            Coverage = coverage;
        }
    }
}