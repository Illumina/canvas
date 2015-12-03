using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace SequencingFiles
{
    /// <summary>
    /// The ReferenceInterval class represents a region of a reference sequence using one-based inclusive start and end positions 
    /// </summary>
    public class ReferenceInterval : IComparable, IComparable<ReferenceInterval>, IEquatable<ReferenceInterval>
    {
        public uint OneBasedStart { get; }
        public uint OneBasedEnd { get; }

        public ReferenceInterval(uint oneBasedStart, uint oneBasedEnd)
        {
            if (oneBasedStart > oneBasedEnd)
                throw new ArgumentException($"Invalid interval. Start position {oneBasedStart} cannot be greater than end position {oneBasedEnd}");
            OneBasedStart = oneBasedStart;
            OneBasedEnd = oneBasedEnd;
        }

        public ReferenceInterval(int oneBasedStart, int oneBasedEnd) : this(checked((uint)oneBasedStart), checked((uint)oneBasedEnd))
        {
        }

        public bool Overlaps(ReferenceInterval interval)
        {
            return !DoesNotOverlap(interval);
        }

        public bool DoesNotOverlap(ReferenceInterval interval)
        {
            return OneBasedEnd < interval.OneBasedStart || interval.OneBasedEnd < OneBasedStart;
        }

        public static List<ReferenceInterval> MergeIntervals(List<ReferenceInterval> intervals)
        {
            List<ReferenceInterval> mergedIntervals = new List<ReferenceInterval>();
            if (!intervals.Any()) return mergedIntervals;

            var sortedIntervals = intervals.OrderBy(r => r);
            var currentInterval = sortedIntervals.First();
            foreach (var nextInterval in sortedIntervals.Skip(1))
            {
                if (currentInterval.Overlaps(nextInterval))
                {
                    uint newEnd = Math.Max(currentInterval.OneBasedEnd, nextInterval.OneBasedEnd);
                    currentInterval = new ReferenceInterval(currentInterval.OneBasedStart, newEnd);
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

        public bool Equals(ReferenceInterval other)
        {
            return
                OneBasedStart == other.OneBasedStart &&
                OneBasedEnd == other.OneBasedEnd;
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
            hash = hash * 31 + OneBasedStart.GetHashCode();
            return hash * 31 + OneBasedEnd.GetHashCode();
        }

        #endregion

        #region IComparable

        int IComparable.CompareTo(object obj)
        {
            if (obj == null)
                throw new ArgumentException("Cannot compare with null");
            var other = obj as ReferenceInterval;
            if (other == null)
                throw new ArgumentException($"{nameof(obj)} must be of the same type");
            return CompareTo(other);
        }

        int IComparable<ReferenceInterval>.CompareTo(ReferenceInterval other)
        {
            return CompareTo(other);
        }

        private int CompareTo(ReferenceInterval other)
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