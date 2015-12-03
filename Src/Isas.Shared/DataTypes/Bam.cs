using Isas.Shared;
using SequencingFiles;
using System;

namespace Illumina.SecondaryAnalysis
{
    public class Bam : IMoveable<Bam>, IMoveableResult<Bam>
    {
        public readonly IFileLocation BamFile;
        private Tuple<int, int?> _readLength;
        private bool? _isPairedEnd;

        public IFileLocation Index { get { return BamFile.AppendName(".bai"); } }

        public IFileLocation MD5 { get { return BamFile.AppendName(".md5"); } }

        public int ReadLength1
        {
            get
            {
                if (_readLength == null)
                    LoadReadLength();
                return _readLength.Item1;
            }
            // No set method makes this property readonly to the outside
        }

        public int? ReadLength2
        {
            get
            {
                if (_readLength == null)
                    LoadReadLength();
                return _readLength.Item2;
            }
            // No set method makes this property readonly to the outside
        }

        public bool IsPairedEnd
        {
            get
            {
                if (_isPairedEnd == null)
                    LoadIsPairedEnd();
                return _isPairedEnd.Value;
            }
            // No set method makes this property readonly to the outside
        }

        // used by json.net to bypass existence check
        private Bam() { }

        public Bam(IFileLocation path)
            : this(path, null, null)
        {
        }

        public Bam(IFileLocation path, int readLength1, int? readLength2 = null)
            : this(path, new Tuple<int, int?>(readLength1, readLength2), readLength2.HasValue)
        {
        }

        private Bam(IFileLocation path, Tuple<int, int?> readLength, bool? isPairedEnd)
        {
            //todo: enable this check as long as no one is depending on using a Bam object to define where a new Bam file should get created
            //if (!path.Exists) throw new ArgumentException(string.Format("Bam file at {0} does not exist", path));
            BamFile = path;
            _readLength = readLength;
            _isPairedEnd = isPairedEnd;
        }

        private void LoadReadLength()
        {
            bool isPairedEnd = false;
            int readLength1 = 0;
            int readLength2 = -1;
            using (BamReader reader = new BamReader(BamFile.FullName))
            {
                BamAlignment read = new BamAlignment();
                for (int readIndex = 0; readIndex < 10000; readIndex++)
                {
                    bool result = reader.GetNextAlignment(ref read, false);
                    if (!result) break;
                    bool currenReadPairedEnd = read.IsPaired();
                    if (currenReadPairedEnd)
                    {
                        isPairedEnd = true;
                    }
                    if (!currenReadPairedEnd || read.IsFirstMate())
                    {
                        readLength1 = Math.Max(read.Bases.Length, readLength1);
                    }
                    else
                    {
                        readLength2 = Math.Max(read.Bases.Length, readLength2);
                    }
                }
            }
            int? readLength2OrNull = null;
            if (readLength2 != -1) readLength2OrNull = readLength2;
            _readLength = new Tuple<int, int?>(readLength1, readLength2OrNull);
            _isPairedEnd = isPairedEnd;
        }

        private void LoadIsPairedEnd()
        {
            _isPairedEnd = false;
            using (BamReader reader = new BamReader(BamFile.FullName))
            {
                BamAlignment read = new BamAlignment();
                if (reader.GetNextAlignment(ref read, false))
                {
                    _isPairedEnd = read.IsPaired();
                }
            }
        }

        public Bam Move(FileNamingConvention newName)
        {
            IFileLocation newBam = BamFile.Move(newName);
            if (Index.Exists)
                Index.Move(newName);
            if (MD5.Exists)
                MD5.Move(newName);
            return new Bam(newBam, _readLength, _isPairedEnd);
        }

        public Bam Move(IFileLocation newPath)
        {
            Utilities.Move(BamFile.FullName, newPath.FullName);
            Bam movedBam = new Bam(newPath, _readLength, _isPairedEnd);
            if (Index.Exists)
                Index.MoveTo(movedBam.Index);
            if (MD5.Exists)
                MD5.MoveTo(movedBam.MD5);
            return movedBam;
        }

        public Bam Move(IFileLocation newPath, Action<IFileLocation, IFileLocation> move)
        {
            move(BamFile, newPath);
            Bam movedBam = new Bam(newPath, _readLength, _isPairedEnd);
            if (Index.Exists)
                move(Index, movedBam.Index);
            if (MD5.Exists)
                move(MD5, movedBam.MD5);
            return movedBam;
        }

        public Bam MoveWithStub(IFileLocation newStub, Action<IFileLocation, IFileLocation> move)
        {
            IFileLocation newPath = newStub.AppendName(".bam");
            return Move(newPath, move);
        }

        public static Bam GetBamFromStub(IFileLocation stub)
        {
            return new Bam(stub.AppendName(".bam"));
        }

        public void Delete()
        {
            BamFile.Delete();
            Index.Delete();
            MD5.Delete();
        }

        public Bam Move(Bam destination)
        {
            BamFile.MoveAndLink(destination.BamFile);
            if (Index.Exists)
                Index.MoveAndLink(destination.Index);
            else if (destination.Index.Exists) // Make sure that if there was a destination index, it is erased
                destination.Index.Delete();
            if (MD5.Exists)
                MD5.MoveAndLink(destination.MD5);
            else if (destination.MD5.Exists) // Make sure that if there was a destination md5, it is erased
                destination.MD5.Delete();
            return destination;
        }
    }

    public static class BamMove
    {
        public static SampleSet<Bam> Move(this SampleSet<Bam> bams, Func<SampleInfo, IFileLocation> getBamPath)
        {
            var newBams = new SampleSet<Bam>();
            foreach (var bam in bams)
            {
                IFileLocation newBamPath = getBamPath(bam.Key);
                newBams[bam.Key] = bam.Value.Move(newBamPath);
            }
            return newBams;
        }
    }
}
