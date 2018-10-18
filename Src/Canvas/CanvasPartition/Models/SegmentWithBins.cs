using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.Statistics;

namespace CanvasPartition.Models
{
    public class SegmentWithBins
    {
        public readonly int Identifier;
        public List<Bin> Bins
        {
            get { return _bins.OrderBy(b => b.Start).ToList(); }
        }
        public uint Start { get; private set; }
        public uint End { get; private set; }
        private readonly List<Bin> _bins;

        public double MedianCoverage
        {
            get { return _bins.Select(b => b.Coverage).Median(); }
        }

        public SegmentWithBins(int identifier, Bin initialBin)
        {
            Start = uint.MaxValue;
            End = uint.MinValue;

            Identifier = identifier;
            _bins = new List<Bin>();

            if (initialBin != null)
            {
                AddBin(initialBin);
            }
        }

        public void AddBin(Bin bin)
        {
            _bins.Add(bin);
            if (bin.Start < Start)
            {
                Start = bin.Start;
            }

            if (bin.End > End)
            {
                End = bin.End;
            }
        }
    }
}