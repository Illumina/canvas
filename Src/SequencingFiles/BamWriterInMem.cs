using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Runtime.InteropServices;
using System.IO;

namespace SequencingFiles
{
    public class BamWriterInMem : BamWriter
    {
        public class BamWriterHandle
        {
            BamWriterInMem _bamWriter;
            int _threadNumber;

            public BamWriterHandle(BamWriterInMem bamWriter, int threadNumber)
            {
                _bamWriter = bamWriter;
                _threadNumber = threadNumber;
            }

            public void WriteAlignment(ref BamAlignment al)
            {
                _bamWriter.WriteAlignment(ref al, _threadNumber);
            }
        }

        #region members
        private List<List<SerializedBamAlignment>> _bamAlignmentLists;
        private readonly long _maxMemPerThread;
        private List<int> _totalMem;
        private int _numTempFiles;
        private ReaderWriterLockSlim _lock = new ReaderWriterLockSlim();
        #endregion

        public BamWriterInMem(string filename,
                              string samHeader,
                              List<GenomeMetadata.SequenceMetadata> references,
                              int maxMemMB = 500, // Same default as samtools sort
                              int numThreads = 1,
                              int compressionLevel = BamConstants.DefaultCompression)
            : base(filename, samHeader, references, compressionLevel, numThreads)
        {
            _bamAlignmentLists = new List<List<SerializedBamAlignment>>();
            _maxMemPerThread = (Convert.ToInt64(maxMemMB)*1048576) / Convert.ToInt64(numThreads);
            _totalMem = new List<int>();
            _numTempFiles = 0;
            _lock = new ReaderWriterLockSlim();

            for (int i = 0; i < numThreads; ++i)
            {
                _bamAlignmentLists.Add(new List<SerializedBamAlignment>());
                _totalMem.Add(0);
            }
        }

        public List<BamWriterHandle> GenerateGetHandles()
        {
            List<BamWriterHandle> handles = new List<BamWriterHandle>();
            for (int i = 0; i < NumThreads; ++i)
            {
                handles.Add(new BamWriterHandle(this, i));
            }
            return handles;
        }

        public void WriteAlignment(ref BamAlignment al)
        {
            WriteAlignment(ref al, 0);
        }

        private void WriteAlignment(ref BamAlignment al, int bufferNumber)
        {
            SerializedBamAlignment serializedAl =
            new SerializedBamAlignment(SerializeAlignment(ref al),
                                        al.RefID,
                                        al.Position,
                                        al.AlignmentFlag,
                                        al.FragmentLength,
                                        al.MapQuality,
                                        al.MatePosition,
                                        al.MateRefID,
                                        al.Name,
                                        al.IsReverseStrand());

            _lock.EnterReadLock();

            // The 1.6 is an empirical value
            _totalMem[bufferNumber] += Convert.ToInt32(1.6*serializedAl.getApproxSize());
            _bamAlignmentLists[bufferNumber].Add(serializedAl);
            _lock.ExitReadLock();

            if (_totalMem[bufferNumber] > _maxMemPerThread)
            {
                _lock.EnterWriteLock();

                ++_numTempFiles;
                SortMemAndWrite();

                _lock.ExitWriteLock();
            }
        }

        public void SortAndWrite()
        {
            if (_numTempFiles > 0)
            {
                ++_numTempFiles;
            }

            SortMemAndWrite();

            if (_numTempFiles > 0)
            {
                MergeTempFiles();
            }
        }

        private string GetTempFileName(int index)
        {
            return Path.Combine(Path.GetDirectoryName(Filename), (index+1).ToString() + Path.GetFileName(Filename));
        }

        // This method is very specialized. Each list should be the same size, so sorting
        // each on its own thread is very efficient. There is one additional pass
        // to merge the lists, but the writing occurs during this iteration, so we're
        // not waiting to write.
        private void SortMemAndWrite()
        {
            string filename = _numTempFiles > 0 ? GetTempFileName(_numTempFiles-1) : Filename;

            // Sort each of the lists
            BamAlignmentComparer comparer = new BamAlignmentComparer();
            Parallel.ForEach(_bamAlignmentLists, bamAlignmentList =>
                {
                    bamAlignmentList.Sort(comparer);
                }
            );

            List<BamRecordAccessor> bamListContainer = new List<BamRecordAccessor>();
            foreach (List<SerializedBamAlignment> bamAlignmentList in _bamAlignmentLists)
            {
                bamListContainer.Add(new BamListContainer(bamAlignmentList));
            }

            BamWriter bamWriter = _numTempFiles > 0 ? new BamWriter(filename, SamHeader, References, CompressionLevel, NumThreads) : this;
            WriteListsInOrder(bamListContainer, bamWriter);

            if (_numTempFiles > 0)
            {
                bamWriter.Dispose();

                for (int i = 0; i < _totalMem.Count; ++i)
                {
                    _totalMem[i] = 0;
                    _bamAlignmentLists[i].Clear();
                }
            }
        }

        #region helper classes

        [StructLayout(LayoutKind.Sequential, Pack = 1)]
        private struct SerializedBamAlignment
        {
            public SerializedBamAlignment(byte[] serializedData, int refID, int position, uint alignmentFlag, int fragmentLength, uint mapQuality, int matePos, int mateRefID, string name, bool isReverseStrand)
            {
                SerializedBam = serializedData;
                RefID = refID;
                Position = position;

                AlignmentFlag = alignmentFlag;
                FragmentLength = fragmentLength;
                MapQuality = mapQuality;
                MatePosition = matePos;
                MateRefID = mateRefID;
                Name = name;
                IsReverseStrand = isReverseStrand;
            }

            public int getApproxSize()
            {
                return SerializedBam.Length + Name.Length + _numBytesValueFields;
            }

            public byte[] SerializedBam;
            public int RefID;
            public int Position;

            public uint AlignmentFlag;
            public int FragmentLength;
            public uint MapQuality;
            public int MatePosition;
            public int MateRefID;
            public string Name;
            public bool IsReverseStrand;
            private const int _numBytesValueFields = 32;
        }

        private class BamAlignmentComparer : IComparer<SerializedBamAlignment>
        {
            // This comparison operator is copied exactly from samtools.
            public int Compare(SerializedBamAlignment a, SerializedBamAlignment b)
            {
                int cmp = FileOrderCompare(a, b);

                if (cmp != 0)
                {
                    return cmp;
                }

                // They're the same.
                // Test of negative strand flag is not really necessary, because it is tested
                // with cmp if getFlags, but it is left here because that is the way it was done
                // in the past.
                if (a.IsReverseStrand == b.IsReverseStrand)
                {
                    cmp = string.CompareOrdinal(a.Name, b.Name);

                    if (cmp != 0) return cmp;
                    cmp = a.AlignmentFlag.CompareTo(b.AlignmentFlag);
                    if (cmp != 0) return cmp;
                    cmp = a.MapQuality.CompareTo(b.MapQuality);
                    if (cmp != 0) return cmp;
                    cmp = a.MateRefID.CompareTo(b.MateRefID);
                    if (cmp != 0) return cmp;
                    cmp = a.MatePosition.CompareTo(b.MatePosition);
                    if (cmp != 0) return cmp;
                    cmp = a.FragmentLength.CompareTo(b.FragmentLength);
                    return cmp;

                }
                else return (a.IsReverseStrand ? 1 : -1);
            }

            private int FileOrderCompare(SerializedBamAlignment a, SerializedBamAlignment b)
            {
                if (a.RefID == -1)
                {
                    return (b.RefID == -1 ? 0 : 1);
                }
                else if (b.RefID == -1)
                {
                    return -1;
                }
                int cmp = a.RefID - b.RefID;
                if (cmp != 0)
                {
                    return cmp;
                }
                return a.Position - b.Position;
            }
        }

        private abstract class BamRecordAccessor
        {
            public abstract SerializedBamAlignment GetCurrentAlignment();
            public abstract void MoveToNextRecord();
            public abstract bool IsEnd();
            public virtual void Close() { }
        }

        private class BamListContainer : BamRecordAccessor
        {
            List<SerializedBamAlignment> _bamList;
            int _index;
            public BamListContainer(List<SerializedBamAlignment> bamList)
            {
                _index = 0;
                _bamList = bamList;
            }

            public override SerializedBamAlignment GetCurrentAlignment()
            {
                return _bamList[_index];
            }

            public override void MoveToNextRecord()
            {
                ++_index;
            }

            public override bool IsEnd()
            {
                return _index >= _bamList.Count;
            }
        }

        private class BamReaderContainer : BamRecordAccessor
        {
            BamReader _bamReader;
            SerializedBamAlignment _currentAlignment;
            bool _isEnd;

            public BamReaderContainer(string filename)
            {
                _bamReader = new BamReader(filename);
                MoveToNextRecord();
            }

            public override SerializedBamAlignment GetCurrentAlignment()
            {
                return _currentAlignment;
            }

            public override void MoveToNextRecord()
            {
                BamAlignment al = new BamAlignment();
                _isEnd = !_bamReader.GetNextAlignment(ref al, false);

                if (_isEnd)
                {
                    return;
                }

                // It's sad that we need to reserialize this, when we just read it.
                _currentAlignment =
                    new SerializedBamAlignment(BamWriter.SerializeAlignment(ref al),
                                               al.RefID,
                                               al.Position,
                                               al.AlignmentFlag,
                                               al.FragmentLength,
                                               al.MapQuality,
                                               al.MatePosition,
                                               al.MateRefID,
                                               al.Name,
                                               al.IsReverseStrand());
            }

            public override bool IsEnd()
            {
                return _isEnd;
            }

            public override void Close()
            {
                base.Close();
                string filePath = _bamReader.BamPath;
                _bamReader.Dispose();

                // Delete the temporary file.
                File.Delete(filePath);
            }
        }

        #endregion

        void WriteListsInOrder(List<BamRecordAccessor> bamContainers, BamWriter bamWriter)
        {
            List<int> bamAlignmentListIndexes = new List<int>();
            for (int i = 0; i < bamContainers.Count; ++i)
            {
                if (!bamContainers[i].IsEnd())
                {
                    bamAlignmentListIndexes.Add(i);
                }
            }

            if (bamAlignmentListIndexes.Count == 0)
            {
                return;
            }

            BamAlignmentComparer comparer = new BamAlignmentComparer();
            while (bamAlignmentListIndexes.Count > 1)
            {
                // Look through all the lists to find the least element
                int whichList = bamAlignmentListIndexes[0];
                int whichIndexList = 0;
                for (int i = 1; i < bamAlignmentListIndexes.Count; ++i)
                {
                    if (0 > comparer.Compare(bamContainers[bamAlignmentListIndexes[i]].GetCurrentAlignment(),
                                             bamContainers[whichList].GetCurrentAlignment()))
                    {
                        whichList = bamAlignmentListIndexes[i];
                        whichIndexList = i;
                    }
                }

                byte[] align = bamContainers[whichList].GetCurrentAlignment().SerializedBam;
                bamWriter.Write(align, Convert.ToUInt32(align.Length));
                bamContainers[whichList].MoveToNextRecord();

                if (bamContainers[whichList].IsEnd())
                {
                    bamAlignmentListIndexes.RemoveAt(whichIndexList);
                }
            }

            // There is only 1 list left
            int listIndex = bamAlignmentListIndexes[0];
            while (!bamContainers[listIndex].IsEnd())
            {
                byte[] al = bamContainers[listIndex].GetCurrentAlignment().SerializedBam;
                bamWriter.Write(al, Convert.ToUInt32(al.Length));
                bamContainers[listIndex].MoveToNextRecord();
            }
        }

        private void MergeTempFiles()
        {
            List<BamRecordAccessor> bamReaders = new List<BamRecordAccessor>();

            for (int i = 0; i < _numTempFiles; ++i)
            {
                string tmpFileName = GetTempFileName(i);

                // Each of these files has already been sorted individually.
                bamReaders.Add(new BamReaderContainer(tmpFileName));
            }

            WriteListsInOrder(bamReaders, this);

            foreach (BamRecordAccessor bamReader in bamReaders)
            {
                bamReader.Close();
            }
        }
    }
}
