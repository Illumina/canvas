using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace SequencingFiles
{
    internal class BamIndexRegion
    {
        public ulong Begin;
        public ulong End;

        public BamIndexRegion(ulong begin, ulong end)
        {
            Begin = begin;
            End = end;
        }
    }

    /// <summary>
    ///     Contains the index offsets and regions for a reference sequence
    /// </summary>
    internal class BamReferenceIndex
    {
        public List<ulong> OffsetList;
        public Dictionary<uint, List<BamIndexRegion>> RegionsDictionary;

        public BamReferenceIndex()
        {
            RegionsDictionary = new Dictionary<uint, List<BamIndexRegion>>();
            OffsetList = new List<ulong>();
        }
    }

    public class BamIndex
    {
        #region member variables

        private const int BamLidxShift = 14;
        public const uint BamMaxBin = 37450;
        private readonly List<BamReferenceIndex> _index;

        private ulong _beginOffset;
        private ulong _endOffset;
        private bool _hasUnalignedReads;
        private uint _lastBin;
        private ulong _lastOffset;
        private int _lastPosition;
        private int _lastRefID;
        private ulong _numAligned;
        private ulong _numUnaligned;
        private ulong _numUnalignedWithoutCoordinates;
        private uint _saveBin;
        private ulong _saveOffset;
        private int _saveRefID;

        public ulong NumUnalignedWithoutCoordinates { get { return _numUnalignedWithoutCoordinates; } }

        #endregion

        // constructor
        public BamIndex()
        {
            _index = new List<BamReferenceIndex>();
        }

        /// <summary>
        ///     Adds a region belonging to a specific bin and reference sequence in the index
        /// </summary>
        private static void AddBamRegion(ref Dictionary<uint, List<BamIndexRegion>> binMap, uint bin, ulong beg,
                                         ulong end)
        {
            List<BamIndexRegion> binList;
            if (binMap.TryGetValue(bin, out binList))
            {
                binList.Add(new BamIndexRegion(beg, end));
            }
            else
            {
                binList = new List<BamIndexRegion> { new BamIndexRegion(beg, end) };
                binMap[bin] = binList;
            }
        }

        /// <summary>
        ///     Adds an offset to a specific reference sequence in the index
        /// </summary>
        private static void AddOffset(ref List<ulong> offsets, ref BamAlignment al, ulong offset)
        {
            int beg = al.Position >> BamLidxShift;
            int end = (al.GetEndPosition() - 1) >> BamLidxShift;

            // initialize additional entries if needed
            while (offsets.Count < (end + 1)) offsets.Add(0);

            if (beg == end)
            {
                if (offsets[beg] == 0) offsets[beg] = offset;
            }
            else
            {
                for (int i = beg; i <= end; i++)
                {
                    if (offsets[i] == 0) offsets[i] = offset;
                }
            }
        }

        /// <summary>
        ///     retrieves the bin keys that overlap the specified position
        /// </summary>
        private static List<uint> GetBinKeysThatOverlapPos(uint begin, uint end)
        {
            // fix the maximum value of the specified position
            List<uint> binKeys = new List<uint>();
            --end;

            // bin '0' always a valid bin
            binKeys.Add(0);

            // get rest of bins that contain this region
            uint k;
            for (k = 1 + (begin >> 26); k <= 1 + (end >> 26); ++k)
            {
                binKeys.Add(k);
            }
            for (k = 9 + (begin >> 23); k <= 9 + (end >> 23); ++k)
            {
                binKeys.Add(k);
            }
            for (k = 73 + (begin >> 20); k <= 73 + (end >> 20); ++k)
            {
                binKeys.Add(k);
            }
            for (k = 585 + (begin >> 17); k <= 585 + (end >> 17); ++k)
            {
                binKeys.Add(k);
            }
            for (k = 4681 + (begin >> 14); k <= 4681 + (end >> 14); ++k)
            {
                binKeys.Add(k);
            }

            // return number of bins stored
            return binKeys;
        }

        /// <summary>
        ///     Merges regions that share start and end in the same bin
        /// </summary>
        private static void MergeRegionList(ref List<BamIndexRegion> regions)
        {
            // merge regions with overlapping bins in the list
            int prevIndex = 0;

            for (int currentIndex = 1; currentIndex < regions.Count; ++currentIndex)
            {
                if ((regions[prevIndex].End >> 16) == (regions[currentIndex].Begin >> 16))
                {
                    regions[prevIndex].End = regions[currentIndex].End;
                }
                else
                {
                    regions[++prevIndex] = regions[currentIndex];
                }
            }

            // remove the remaining regions
            int numMergedRegions = prevIndex + 1;
            int numRegionsToDelete = regions.Count - numMergedRegions;
            if (numRegionsToDelete > 0) regions.RemoveRange(numMergedRegions, numRegionsToDelete);
        }

        /// <summary>
        ///     Merges the regions that occur in each reference sequence and bin
        /// </summary>
        private void MergeBinnedRegions()
        {
            // iterate over indices for each reference sequence
            foreach (BamReferenceIndex t in _index)
            {
                // grab the keys
                Dictionary<uint, List<BamIndexRegion>> regionDictionary = t.RegionsDictionary;
                uint[] saveBins = new uint[regionDictionary.Keys.Count];
                regionDictionary.Keys.CopyTo(saveBins, 0);

                // iterate over the bins
                foreach (uint t1 in saveBins)
                {
                    List<BamIndexRegion> regions = regionDictionary[t1];
                    MergeRegionList(ref regions);
                    regionDictionary[t1] = regions;
                }
            }
        }

        /// <summary>
        ///     Fills missing offsets with the preceding non-zero offset
        /// </summary>
        private void FillMissing()
        {
            foreach (BamReferenceIndex t in _index)
            {
                List<ulong> offsetList = t.OffsetList;
                for (int i = 1; i < offsetList.Count; i++)
                {
                    if (offsetList[i] == 0) offsetList[i] = offsetList[i - 1];
                }
            }
        }

        /// <summary>
        ///     Creates an index from a specified BAM file
        /// </summary>
        public void CreateIndexFromBamFile(string filename)
        {
            _numUnalignedWithoutCoordinates = 0;

            // open the BAM file and retrieve the reference data
            using (BamReader reader = new BamReader(filename))
            {
                // allocate space for the reference index
                List<GenomeMetadata.SequenceMetadata> references = reader.GetReferences();

                // iterate over all of the reads in the BAM file
                Initialize(references.Count, reader.Tell());

                BamAlignment al = new BamAlignment();
                while (reader.GetNextAlignment(ref al, true))
                {
                    if (!UpdateReferenceIndex(ref al, reader.Tell())) break;
                }

                // perform some post-processing on the index
                PostProcessing(reader.Tell());

                if (_hasUnalignedReads)
                {
                    while (reader.GetNextAlignment(ref al, true)) ++_numUnalignedWithoutCoordinates;
                }
            }

            // write the index to a file
            WriteIndex(filename + ".bai");
        }

        /// <summary>
        ///     returns a list of the index regions for the desire reference sequence and position
        /// </summary>
        internal bool GetOffsets(int refID, int position, out BamIterator bamIterator)
        {
            // initialize the bam iterator
            bamIterator = new BamIterator(refID, position, 1 << 29);

            // adjust the specified position
            if (position < 0) position = 0;

            // calculate which bins overlap this region
            List<uint> binKeys = GetBinKeysThatOverlapPos((uint)bamIterator.Begin, (uint)bamIterator.End);
            // get bins and offsets for this reference
            BamReferenceIndex refIndex = _index[refID];
            Dictionary<uint, List<BamIndexRegion>> binnedRegions = refIndex.RegionsDictionary;
            List<ulong> offsets = refIndex.OffsetList;
            ulong minOffset = 0;

            if (offsets.Count > 0)
            {
                int offsetIndex = position >> BamLidxShift;
                minOffset = (offsetIndex >= offsets.Count) ? offsets[offsets.Count - 1] : offsets[offsetIndex];

                // improvement for index files built by tabix prior to 0.1.4
                if (minOffset == 0)
                {
                    // Scan backward for a valid offset:
                    if (offsetIndex > offsets.Count) offsetIndex = offsets.Count;

                    int i;
                    for (i = offsetIndex - 1; i >= 0; i--) if (offsets[i] != 0) break;
                    if (i >= 0) minOffset = offsets[i];
                }
            }

            // get the total count of regions represented in each bin
            // this is a check to see if we can exit early
            int numIndexRegions = 0;
            List<BamIndexRegion> bamIndexRegions;

            foreach (uint binKey in binKeys)
            {
                if (binnedRegions.TryGetValue(binKey, out bamIndexRegions))
                {
                    numIndexRegions += bamIndexRegions.Count;
                }
            }

            if (numIndexRegions == 0) return false;

            // grab all of the index regions that end after the minimum offset
            List<BamIndexRegion> regionsAfterMinOffset = new List<BamIndexRegion>(numIndexRegions);

            foreach (uint binKey in binKeys)
            {
                if (binnedRegions.TryGetValue(binKey, out bamIndexRegions))
                {
                    foreach (BamIndexRegion indexRegion in bamIndexRegions)
                    {
                        if (indexRegion.End > minOffset) regionsAfterMinOffset.Add(indexRegion);
                    }
                }
            }

            // sort the index regions
            BamIndexRegion[] sortedList = regionsAfterMinOffset.OrderBy(x => x.Begin).ThenBy(x => x.End).ToArray();

            // consolidate completely contained adjacent blocks
            int prevIndex = 0;
            for (int currentIndex = 1; currentIndex < sortedList.Length; ++currentIndex)
            {
                if (sortedList[prevIndex].End < sortedList[currentIndex].End)
                    sortedList[++prevIndex] = sortedList[currentIndex];
            }

            numIndexRegions = prevIndex + 1;

            // resolve overlaps between adjacent blocks; this may happen due to the merge in indexing
            for (int currentIndex = 1; currentIndex < numIndexRegions; ++currentIndex)
            {
                if (sortedList[currentIndex - 1].End >= sortedList[currentIndex].Begin)
                {
                    sortedList[currentIndex - 1].End = sortedList[currentIndex].Begin;
                }
            }

            // merge adjacent blocks
            prevIndex = 0;

            for (int currentIndex = 1; currentIndex < numIndexRegions; ++currentIndex)
            {
                if ((sortedList[prevIndex].End >> 16) == (sortedList[currentIndex].Begin >> 16))
                {
                    sortedList[prevIndex].End = sortedList[currentIndex].End;
                }
                else
                {
                    sortedList[++prevIndex] = sortedList[currentIndex];
                }
            }

            numIndexRegions = prevIndex + 1;

            // add the index regions to our list
            bamIterator.Offsets = new BamIndexRegion[numIndexRegions];
            for (int currentIndex = 0; currentIndex < numIndexRegions; ++currentIndex)
            {
                bamIterator.Offsets[currentIndex] = sortedList[currentIndex];
            }

            return true;
        }

        /// <summary>
        /// returns the largest BAM offset in the specified reference sequence
        /// </summary>
        public ulong GetLargestBamOffset()
        {
            // initialize
            long largestBlockAddress = -1;
            ulong largestBamOffset   = 0;

            // check all of the offsets for this reference sequence
            for (int targetIndex = _index.Count - 1; targetIndex >= 0; --targetIndex)
            {
                foreach (var offset in _index[targetIndex].OffsetList)
                {
                    // calculate the block address
                    long blockAddress = (long)((offset >> 16) & 0xFFFFFFFFFFFF);

                    if (blockAddress > largestBlockAddress)
                    {
                        largestBlockAddress = blockAddress;
                        largestBamOffset    = offset;
                    }
                }

                if (largestBlockAddress != -1) break;
            }

            return largestBamOffset;
        }

        public void Initialize(int numReferences, ulong offset)
        {
            _lastRefID = int.MinValue;
            _saveRefID = int.MinValue;
            _lastBin = uint.MaxValue;
            _saveBin = uint.MaxValue;
            _lastPosition = int.MinValue;

            _saveOffset = offset;
            _lastOffset = offset;
            _beginOffset = offset;
            _endOffset = offset;

            _numAligned = 0;
            _numUnaligned = 0;

            _hasUnalignedReads = false;

            _index.Clear();
            for (int i = 0; i < numReferences; i++)
            {
                _index.Add(new BamReferenceIndex());
            }
        }

        /// <summary>
        ///     Performs some final processing on the index before it can be written to a file
        /// </summary>
        public void PostProcessing(ulong offset)
        {
            if (_saveRefID >= 0)
            {
                AddBamRegion(ref _index[_saveRefID].RegionsDictionary, _saveBin, _saveOffset, offset);
                AddBamRegion(ref _index[_saveRefID].RegionsDictionary, BamMaxBin, _beginOffset, offset);
                AddBamRegion(ref _index[_saveRefID].RegionsDictionary, BamMaxBin, _numAligned, _numUnaligned);
            }

            MergeBinnedRegions();
            FillMissing();
        }

        /// <summary>
        ///     Reads a BAM index from the supplied filename
        /// </summary>
        /// <returns>true if the index was successfully loaded</returns>
        public bool ReadIndex(string filename)
        {
            // check if the file exists
            if (!File.Exists(filename)) return false;

            using (
                BinaryReader reader =
                    new BinaryReader(new FileStream(filename, FileMode.Open, FileAccess.Read, FileShare.Read)))
            {
                // check to see if we have a proper BAM index signature
                byte[] buffer = reader.ReadBytes(4);

                string magicNumberString = Encoding.ASCII.GetString(buffer, 0, 4);

                if (magicNumberString != BamConstants.BaiMagicNumber)
                {
                    throw new ApplicationException(
                        string.Format("ERROR: Expected the BAM index magic number to be {0}, but found {1}.",
                                      BamConstants.BaiMagicNumber, magicNumberString));
                }

                // read the number of reference sequences
                uint numReferenceSequences = reader.ReadUInt32();

                // iterate over each reference sequence
                _index.Clear();

                for (uint refSeqIndex = 0; refSeqIndex < numReferenceSequences; ++refSeqIndex)
                {
                    // =======================
                    // read the binning index
                    // =======================

                    BamReferenceIndex refIndex = new BamReferenceIndex();

                    // read the number of bins
                    uint numBins = reader.ReadUInt32();

                    // iterate over each bin in the regions dictionary
                    for (uint binIndex = 0; binIndex < numBins; ++binIndex)
                    {
                        // read the bin and the number of index regions
                        uint bin = reader.ReadUInt32();
                        uint numIndexRegions = reader.ReadUInt32();

                        // read all of the index regions
                        List<BamIndexRegion> bamIndexRegions = new List<BamIndexRegion>();

                        for (uint regionIndex = 0; regionIndex < numIndexRegions; ++regionIndex)
                        {
                            ulong begin = reader.ReadUInt64();
                            ulong end = reader.ReadUInt64();
                            BamIndexRegion indexRegion = new BamIndexRegion(begin, end);
                            bamIndexRegions.Add(indexRegion);
                        }

                        // add the index regions to our dictionary
                        refIndex.RegionsDictionary[bin] = bamIndexRegions;
                    }

                    // =====================
                    // read the linear index
                    // =====================

                    // read the linear index size
                    uint numOffsets = reader.ReadUInt32();

                    for (uint offsetIndex = 0; offsetIndex < numOffsets; ++offsetIndex)
                    {
                        refIndex.OffsetList.Add(reader.ReadUInt64());
                    }

                    // add the reference index
                    _index.Add(refIndex);
                }

                // read the number of unaligned reads without coordinates
                _numUnalignedWithoutCoordinates = reader.ReadUInt64();
            }

            return true;
        }

        /// <summary>
        ///     Updates the index with respect to the current alignment
        /// </summary>
        /// <returns>false if multiple reads without coordinates are encountered</returns>
        public bool UpdateReferenceIndex(ref BamAlignment alignment, ulong offset)
        {
            // record the number of unaligned reads
            if (alignment.RefID < 0) ++_numUnalignedWithoutCoordinates;

            // update the reference IDs and check that the alignment is sorted
            if (alignment.RefID != _lastRefID)
            {
                _lastRefID = alignment.RefID;
                _lastBin = int.MaxValue;
            }
            else if (alignment.Position < _lastPosition)
            {
                throw new ApplicationException(
                    string.Format(
                        "ERROR: The BAM file is not sorted. An alignment ({0}) occurred before the preceding alignment ({1}).",
                        alignment.Position, _lastPosition));
            }

            if (alignment.RefID >= 0) AddOffset(ref _index[alignment.RefID].OffsetList, ref alignment, _lastOffset);

            if (alignment.Bin != _lastBin)
            {
                if (_saveBin != uint.MaxValue)
                    AddBamRegion(ref _index[_saveRefID].RegionsDictionary, _saveBin, _saveOffset, _lastOffset);
                if ((_lastBin == uint.MaxValue) && (_saveRefID != int.MinValue))
                {
                    _endOffset = _lastOffset;
                    AddBamRegion(ref _index[_saveRefID].RegionsDictionary, BamMaxBin, _beginOffset, _endOffset);
                    AddBamRegion(ref _index[_saveRefID].RegionsDictionary, BamMaxBin, _numAligned, _numUnaligned);
                    _numAligned = _numUnaligned = 0;
                    _beginOffset = _endOffset;
                }

                _saveOffset = _lastOffset;
                _saveBin = _lastBin = alignment.Bin;
                _saveRefID = alignment.RefID;

                if (_saveRefID < 0)
                {
                    _hasUnalignedReads = true;
                    return false;
                }
            }

            if (offset <= _lastOffset)
            {
                throw new ApplicationException(
                    "ERROR: While updating the BAM index, the offset did not increase after processing the last alignment.");
            }

            if (alignment.IsMapped()) ++_numAligned;
            else ++_numUnaligned;

            _lastOffset = offset;
            _lastPosition = alignment.Position;

            return true;
        }

        public void WriteIndex(string filename)
        {
            using (BinaryWriter writer = new BinaryWriter(new FileStream(filename, FileMode.Create)))
            {
                // write the BAM index signature
                writer.Write(Encoding.ASCII.GetBytes(BamConstants.BaiMagicNumber));

                // write the number of reference sequences
                uint numReferenceSequences = (uint)_index.Count;
                writer.Write(BitConverter.GetBytes(numReferenceSequences));

                // iterate over each reference sequence
                foreach (BamReferenceIndex refIndex in _index)
                {
                    // =======================
                    // write the binning index
                    // =======================

                    // write the binning index size
                    uint numBins = (uint)refIndex.RegionsDictionary.Count;
                    writer.Write(BitConverter.GetBytes(numBins));

                    // iterate over each bin in the regions dictionary
                    foreach (KeyValuePair<uint, List<BamIndexRegion>> kvp in refIndex.RegionsDictionary)
                    {
                        // write the bin
                        writer.Write(BitConverter.GetBytes(kvp.Key));

                        // write the number of index regions
                        uint numIndexRegions = (uint)kvp.Value.Count;
                        writer.Write(BitConverter.GetBytes(numIndexRegions));

                        // write out all regions
                        foreach (BamIndexRegion region in kvp.Value)
                        {
                            writer.Write(BitConverter.GetBytes(region.Begin));
                            writer.Write(BitConverter.GetBytes(region.End));
                        }
                    }

                    // ======================
                    // write the linear index
                    // ======================

                    // write the linear index size
                    uint numOffsets = (uint)refIndex.OffsetList.Count;
                    writer.Write(BitConverter.GetBytes(numOffsets));

                    foreach (ulong offset in refIndex.OffsetList)
                    {
                        writer.Write(BitConverter.GetBytes(offset));
                    }
                }

                // write the number of unaligned reads without coordinates
                writer.Write(BitConverter.GetBytes(_numUnalignedWithoutCoordinates));
            }
        }

        /// <summary>
        /// returns a string representation of the BAM index. This produces output
        /// similar to that of idxstats
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            for (int refIndex = 0; refIndex < _index.Count; refIndex++)
            {
                var currentRefIndex = _index[refIndex];
                if (currentRefIndex.RegionsDictionary.Count == 0)
                {
                    sb.AppendFormat("{0}\t0\t0\n", refIndex);
                    continue;
                }

                List<BamIndexRegion> bamIndexRegions;
                if (currentRefIndex.RegionsDictionary.TryGetValue(BamMaxBin, out bamIndexRegions))
                {
                    sb.AppendFormat("{0}\t{1}\t{2}\n", refIndex, bamIndexRegions[1].Begin, bamIndexRegions[1].End);
                }
                else
                {
                    sb.AppendFormat("{0}\t0\t0\n", refIndex);
                }
            }

            sb.AppendFormat("*\t0\t{0}\n", _numUnalignedWithoutCoordinates);

            return sb.ToString();
        }
    }

    internal sealed class BamIterator
    {
        public int Begin;
        public ulong CurrentOffset;
        public int End;
        public BamIndexRegion[] Offsets;
        public int RefID;

        // constructor
        public BamIterator(int refID, int begin, int end)
        {
            RefID = refID;
            Begin = begin;
            End = end;
            CurrentOffset = 0;
            Offsets = null;
        }
    }
}