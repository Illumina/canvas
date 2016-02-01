using System;
using System.Collections.Generic;
using ProtoBuf;

namespace Isas.Shared
{
    [Serializable]
    [ProtoContract]
    public class GenomicInterval : IComparable
    {
        [ProtoMember(1)]
        public int End;
        [ProtoMember(2)]
        public int Start;

        // Parameterless constructor, to make protobuf happy:
        public GenomicInterval() { }

        public GenomicInterval(int start, int end)
        {
            Start = start;
            End = end;
        }

        int IComparable.CompareTo(object obj)
        {
            GenomicInterval otherInterval = (GenomicInterval)obj;
            if (Start < otherInterval.Start) return -1;
            if (Start > otherInterval.Start) return 1;
            return 0;
        }

        /// <summary>
        ///     Return a sorted list of intervals, merging overlapping intervals.
        /// </summary>
        public static void MergeIntervalList(List<GenomicInterval> intervals)
        {
            intervals.Sort();
            int index = 0;
            while (index < intervals.Count - 1)
            {
                if (intervals[index].End > intervals[index + 1].Start)
                {
                    intervals[index].End = Math.Max(intervals[index].End, intervals[index + 1].End);
                    intervals.RemoveAt(index + 1);
                }
                else
                {
                    index++;
                }
            }
        }
    }
}