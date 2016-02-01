using System;
using System.Collections.Generic;
using System.Linq;
using ILMNcommon.Common;

namespace Illumina.SecondaryAnalysis
{
    public class SampleIndex
    {
        public string Index1Sequence { get; }
        public string Index2Sequence { get; }

        public SampleIndex(string index1Sequence, string index2Sequence)
        {
            // treat null same as empty
            Index1Sequence = index1Sequence ?? "";
            Index2Sequence = index2Sequence ?? "";
        }

        public static List<SampleIndex> AllIndexes(SampleSet<HashSet<SampleIndex>> laneIndexes)
        {
            return laneIndexes.SampleData.SelectMany(index => index).ToList();
        }

        private static int MismatchCount(string indexA, string indexB)
        {
            if (indexA.IsNullOrEmpty() || indexB.IsNullOrEmpty()) return 0;

            int mismatchCount = 0;
            for (int position = 0; position < indexA.Length; position++)
            {
                char a = indexA[position];
                char b = indexB[position];
                if (a != 'N' && b != 'N' && a != b) mismatchCount++;
            }
            return mismatchCount;
        }

        public int MismatchCount(SampleIndex otherIndex)
        {
            return
                MismatchCount(Index1Sequence, otherIndex.Index1Sequence) +
                MismatchCount(Index2Sequence, otherIndex.Index2Sequence);
        }

        public SampleIndex Truncate(int index1Length, int index2Length)
        {
            var truncatedIndex1 = Index1Sequence.Length > index1Length ? Index1Sequence.Substring(0, index1Length) : Index1Sequence;
            var truncatedIndex2 = Index2Sequence.Length > index2Length ? Index2Sequence.Substring(0, index2Length) : Index2Sequence;
            return new SampleIndex(truncatedIndex1, truncatedIndex2);
        }

        public bool Empty()
        {
            return Index1Sequence.IsNullOrEmpty() && Index2Sequence.IsNullOrEmpty();
        }

        public string ToString(string indexSeparator)
        {
            string result = "";
            if (!Index1Sequence.IsNullOrEmpty())
                result += Index1Sequence;
            if (!Index2Sequence.IsNullOrEmpty())
            {
                if (!Index1Sequence.IsNullOrEmpty())
                    result += indexSeparator;
                result += Index2Sequence;
            }
            return result;
        }

        /// <summary>
        /// Truncate all indexes at once so we output a unique set of indexes after truncating 
        /// </summary>
        /// <param name="indexes"></param>
        /// <param name="minIndex1Length"></param>
        /// <param name="minIndex2Length"></param>
        /// <returns></returns>
        public static HashSet<SampleIndex> Truncate(HashSet<SampleIndex> indexes, int minIndex1Length, int minIndex2Length)
        {
            HashSet<SampleIndex> truncatedIndexes = new HashSet<SampleIndex>();
            foreach (var index in indexes)
            {
                var truncatedIndex = index.Truncate(minIndex1Length, minIndex2Length);
                truncatedIndexes.Add(truncatedIndex);
            }
            return truncatedIndexes;
        }

        public string LengthString(string indexLengthSeparator = ",")
        {
            string result = "";
            if (!Index1Sequence.IsNullOrEmpty())
                result += Index1Sequence.Length;
            if (!Index2Sequence.IsNullOrEmpty())
            {
                if (!Index1Sequence.IsNullOrEmpty())
                    result += indexLengthSeparator;
                result += Index2Sequence.Length;
            }
            if (result.IsNullOrEmpty()) return "0";
            return result;
        }

        public int TotalLength()
        {
            int result = 0;
            if (!Index1Sequence.IsNullOrEmpty())
                result += Index1Sequence.Length;
            if (!Index2Sequence.IsNullOrEmpty())
            {
                result += Index2Sequence.Length;
            }
            return result;
        }

        public override string ToString()
        {
            return ToString(",");
        }

        public static IEnumerable<Tuple<KeyValuePair<SampleInfo, SampleIndex>, KeyValuePair<SampleInfo, SampleIndex>>> GetIndexCombinations(SampleSet<HashSet<SampleIndex>> sampleIndexes)
        {
            foreach (var indexPair in sampleIndexes.GetPairs())
            {
                foreach (var indexItem1 in indexPair.Item1.Value)
                {
                    foreach (var indexItem2 in indexPair.Item2.Value)
                    {
                        yield return Tuple.Create(new KeyValuePair<SampleInfo, SampleIndex>(indexPair.Item1.Key, indexItem1), new KeyValuePair<SampleInfo, SampleIndex>(indexPair.Item2.Key, indexItem2));
                    }
                }
            }
        }

        public override bool Equals(System.Object obj)
        {
            // If parameter is null return false.
            if (obj == null)
            {
                return false;
            }

            SampleIndex p = obj as SampleIndex;
            if (p == null)
            {
                return false;
            }

            // Return true if the fields match:
            return Index1Sequence == p.Index1Sequence && Index2Sequence == p.Index2Sequence;
        }

        public override int GetHashCode()
        {
            // index1
            return (Index1Sequence + Index2Sequence).GetHashCode();
        }

        public SampleIndex Index1Only()
        {
            return new SampleIndex(Index1Sequence, null);
        }

        public SampleIndex Index2Only()
        {
            return new SampleIndex(null, Index2Sequence);
        }
    }
}