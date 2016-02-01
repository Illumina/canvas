using System;
using System.Collections.Generic;
using System.Linq;

namespace SequencingFiles
{
    /// <summary>
    /// The Interval class represents a span of sequence positions using one-based inclusive start and end positions 
    /// </summary>
    public class Interval : IComparable, IComparable<Interval>, IEquatable<Interval>
    {
        public int OneBasedStart { get; }
        public int OneBasedEnd { get; }

        public Interval(uint oneBasedStart, uint oneBasedEnd) : this(unchecked((int)oneBasedStart), unchecked((int)oneBasedEnd))
        {
        }

        public Interval Shift(int positions)
        {
            return new Interval(OneBasedStart + positions, OneBasedEnd + positions);
        }

        public Interval(int oneBasedStart, int oneBasedEnd)
        {
            if (oneBasedStart < 0 || oneBasedEnd < 0)
                throw new ArgumentException("Invalid interval. Start and end must be greater than zero");
            if (oneBasedStart > oneBasedEnd)
                throw new ArgumentException($"Invalid interval. Start position {oneBasedStart} cannot be greater than end position {oneBasedEnd}");

            OneBasedStart = oneBasedStart;
            OneBasedEnd = oneBasedEnd;
        }

        public bool Overlaps(Interval interval)
        {
            return !DoesNotOverlap(interval);
        }

        public bool DoesNotOverlap(Interval interval)
        {
            return OneBasedEnd < interval.OneBasedStart || interval.OneBasedEnd < OneBasedStart;
        }

        public static List<Interval> MergeIntervals(List<Interval> intervals)
        {
            List<Interval> mergedIntervals = new List<Interval>();
            if (!intervals.Any()) return mergedIntervals;

            var sortedIntervals = intervals.OrderBy(r => r);
            var currentInterval = sortedIntervals.First();
            foreach (var nextInterval in sortedIntervals.Skip(1))
            {
                if (currentInterval.Overlaps(nextInterval))
                {
                    int newEnd = Math.Max(currentInterval.OneBasedEnd, nextInterval.OneBasedEnd);
                    currentInterval = new Interval(currentInterval.OneBasedStart, newEnd);
                }
                else
                {
                    mergedIntervals.Add(currentInterval);
                    currentInterval = nextInterval;
                }
            }
            mergedIntervals.Add(currentInterval);
            return mergedIntervals;
        }

        #region IEquatable

        public bool Equals(Interval other)
        {
            return
                OneBasedStart == other.OneBasedStart &&
                OneBasedEnd == other.OneBasedEnd;
        }

        public override bool Equals(object obj)
        {
            var item = obj as Interval;
            if (item == null) return false;
            return Equals(item);
        }

        public override int GetHashCode()
        {
            int hash = 23;
            hash = hash * 31 + OneBasedStart.GetHashCode();
            return hash * 31 + OneBasedEnd.GetHashCode();
        }

        #endregion

        #region IComparable

        public int CompareTo(object obj)
        {
            if (obj == null)
                throw new ArgumentException("Cannot compare with null");
            var other = obj as Interval;
            if (other == null)
                throw new ArgumentException($"{nameof(obj)} must be of the same type");
            return CompareTo(other);
        }

        public int CompareTo(Interval other)
        {
            if (other == null)
                throw new ArgumentException("Cannot compare with null");
            if (OneBasedStart < other.OneBasedStart) return -1;
            if (OneBasedStart > other.OneBasedStart) return 1;
            if (OneBasedEnd < other.OneBasedEnd) return -1;
            if (OneBasedEnd > other.OneBasedEnd) return 1;
            return 0;
        }

        #endregion
    }
}