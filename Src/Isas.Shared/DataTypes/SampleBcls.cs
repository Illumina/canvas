using System;
using System.Collections.Generic;
using System.Linq;
using Isas.Shared;
using Isas.Shared.Utilities;

namespace Illumina.SecondaryAnalysis
{
    public class SampleBcls
    {
        public readonly BclRunFolder BclRunFolder;
        public readonly Dictionary<int, HashSet<SampleIndex>> LaneIndexes; //only demux/convert reads in these lanes

        private SampleBcls(BclRunFolder bclRunFolder, Dictionary<int, HashSet<SampleIndex>> laneIndexes)
        {
            LaneIndexes = laneIndexes;
            BclRunFolder = bclRunFolder;
            foreach (int lane in LaneIndexes.Keys)
            {
                foreach (SampleIndex index in LaneIndexes[lane])
                {
                    string index1 = index.Index1Sequence;
                    string index2 = index.Index2Sequence;
                    if (BclRunFolder.ReadStructure.Index1 == null && !string.IsNullOrEmpty(index1))
                        throw new ArgumentException("Cannot have an index 1 sequence when the BCL RunFolder does not have an index 1 read");
                    if (BclRunFolder.ReadStructure.Index1 != null && !string.IsNullOrEmpty(index1) && index1.Length > BclRunFolder.ReadStructure.Index1.Length)
                        throw new ArgumentException(
                            $"Sample's index 1 sequence length ({index1.Length}) cannot be greater than the BCL RunFolder index 1 read length ({BclRunFolder.ReadStructure.Index1.Length})");
                    if (BclRunFolder.ReadStructure.Index2 == null && !string.IsNullOrEmpty(index2))
                        throw new ArgumentException("Cannot have an index 2 sequence when the BCL RunFolder does not have an index 2 read");
                    if (BclRunFolder.ReadStructure.Index2 != null && !string.IsNullOrEmpty(index2) && index2.Length > BclRunFolder.ReadStructure.Index2.Length)
                        throw new ArgumentException(
                            $"Sample's index 2 sequence length ({index2.Length}) cannot be greater than the BCL RunFolder index 2 read length ({BclRunFolder.ReadStructure.Index2.Length})");
                }
            }
        }

        public class Builder
        {
            private readonly BclRunFolder _runFolder;
            private readonly Dictionary<int, HashSet<SampleIndex>> _laneIndexes = new Dictionary<int, HashSet<SampleIndex>>();
            public Builder(BclRunFolder runFolder)
            {
                _runFolder = runFolder;
            }

            public void AddIndex(int lane, SampleIndex index)
            {
                HashSet<SampleIndex> laneIndexes;
                if (!_laneIndexes.TryGetValue(lane, out laneIndexes))
                {
                    laneIndexes = new HashSet<SampleIndex>();
                    _laneIndexes[lane] = laneIndexes;
                }
                _laneIndexes[lane].Add(index);
            }

            public bool ContainsLane(int lane)
            {
                return _laneIndexes.ContainsKey(lane);
            }

            public bool ContainsIndex(int lane, SampleIndex index)
            {
                return ContainsLane(lane) && _laneIndexes[lane].Contains(index);
            }

            public IEnumerable<SampleIndex> Indexes(int lane)
            {
                return ContainsLane(lane) ? _laneIndexes[lane] : Enumerable.Empty<SampleIndex>();
            }

            public SampleBcls Create()
            {
                return new SampleBcls(_runFolder, _laneIndexes);
            }
        }

        public static int MinSampleIndex1Length(SampleSet<SampleBcls> bcls)
        {
            var sampleIndexes = SampleIndexes(bcls);
            return MinSampleIndex1Length(sampleIndexes);
        }

        private static IEnumerable<SampleIndex> SampleIndexes(SampleSet<SampleBcls> bcls)
        {
            return bcls.SampleData.SelectMany(sampleBcls => sampleBcls.LaneIndexes.Values.SelectMany(indexes => indexes));
        }

        public static int MinSampleIndex1Length(IEnumerable<SampleIndex> indexes)
        {
            return indexes.Select(index => index.Index1Sequence).MinLength();
        }

        public static int MinSampleIndex2Length(SampleSet<SampleBcls> bcls)
        {
            var sampleIndexes = SampleIndexes(bcls);
            return MinSampleIndex2Length(sampleIndexes);
        }

        public static int MinSampleIndex2Length(IEnumerable<SampleIndex> indexes)
        {
            return indexes.Select(index => index.Index2Sequence).MinLength();
        }
        public static SampleSet<SampleBcls> GetSampleBclsForRunFolder(SampleSet<IEnumerable<SampleBcls>> sampleBcls, BclRunFolder bclRunFolder)
        {
            return sampleBcls.SelectData(bcls => bcls.First(
                b => b.BclRunFolder.RunFolder.FullName.Equals(bclRunFolder.RunFolder.FullName, Utilities.IsThisMono() ? StringComparison.Ordinal : StringComparison.OrdinalIgnoreCase)));
        }

        public static IEnumerable<IDirectoryLocation> GetRunFolders(SampleSet<IEnumerable<SampleBcls>> samples)
        {
            return GetBclRunFolders(samples).Select(b => b.RunFolder);
        }

        public static List<BclRunFolder> GetBclRunFolders(SampleSet<IEnumerable<SampleBcls>> inputSamples)
        {
            return inputSamples.SampleData.SelectMany(bcls => bcls.Select(b => b.BclRunFolder))
                .Distinct().ToList();
        }

        public static IEnumerable<int> GetLanesUsed(SampleSet<SampleBcls> bcls)
        {
            return bcls.SampleData.SelectMany(sampleBcls => sampleBcls.LaneIndexes.Keys).Distinct();
        }

        public static SampleSet<HashSet<SampleIndex>> GetSampleIndexesForLane(SampleSet<SampleBcls> bcls, int lane)
        {
            return bcls
                .WhereData(sampleBcls => sampleBcls.LaneIndexes.ContainsKey(lane))
                .SelectData(sampleBcls => sampleBcls.LaneIndexes[lane]);
        }
    }
}
