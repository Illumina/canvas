using System;

namespace SequencingFiles
{
    /// <summary>
    /// The ReferenceInterval class represents a region of a reference sequence using one-based inclusive start and end positions 
    /// </summary>
    public class ReferenceInterval : IEquatable<ReferenceInterval>
    {
        public string Chromosome { get; }
        public Interval Interval { get; }

        public ReferenceInterval(string chromosome, Interval interval)
        {
            Chromosome = chromosome;
            Interval = interval;
        }

        public ReferenceInterval(ReferencePosition start, ReferencePosition end)
        {
            if (start.Chromosome != end.Chromosome)
                throw new ArgumentException("Start and end positions must be on the same chromosome.");
            Chromosome = start.Chromosome;
            Interval = new Interval(start.Position, end.Position);
        }

        public bool Overlaps(ReferenceInterval other)
        {
            if (Chromosome != other.Chromosome) return false;
            return Interval.Overlaps(other.Interval);
        }

        public override string ToString()
        {
            return $"{Chromosome}:{Interval.OneBasedStart}-{Interval.OneBasedEnd}";
        }

        #region IEquatable

        public bool Equals(ReferenceInterval other)
        {
            return
                other.Chromosome == Chromosome &&
                other.Interval.Equals(Interval);
        }

        public override bool Equals(object obj)
        {
            var item = obj as ReferenceInterval;
            if (item == null) return false;
            return Equals(item);
        }

        public override int GetHashCode()
        {
            int hash = 23;
            hash = hash * 31 + Chromosome.GetHashCode();
            return hash * 31 + Interval.GetHashCode();
        }

        #endregion
    }

    /// <summary>
    /// The ReferencePosition class represents a single one-based position in a reference sequence 
    /// </summary>
    public class ReferencePosition : ReferenceInterval
    {
        public int Position => Interval.OneBasedStart;

        public ReferencePosition(string chromosome, int position) : base(chromosome, new Interval(position, position))
        {
        }

        public ReferencePosition Next()
        {
            return Shift(1);
        }

        public ReferencePosition Shift(int positions)
        {
            return new ReferencePosition(Chromosome, Position + positions);
        }

        /// <summary>
        /// Return the interval spanned by this position and the position after shifting the specified number of positions
        /// </summary>
        /// <param name="shiftedPositions">
        /// the number of positions to shift this position. This defines the endpoint for the returned interval. 
        /// A negative value is valid provided that the starting position of the resulting interval is not negative
        /// </param>
        /// <returns>
        /// the interval spanned by this position and the position after shifting the specified number of positions. 
        /// the length of the resulting interval will be one greater than the specified number of positions
        /// </returns>
        public ReferenceInterval Span(int shiftedPositions)
        {
            return shiftedPositions > 0 ? new ReferenceInterval(this, Shift(shiftedPositions)) : new ReferenceInterval(Shift(shiftedPositions), this);
        }

        public override string ToString()
        {
            return $"{Chromosome}:{Position}";
        }

        public ReferencePosition Previous()
        {
            return Shift(-1);
        }
    }
}