using System;
using System.IO;
using Illumina.Common.FileSystem;

namespace CanvasTest
{
    /// <summary>
    /// A fixture for testing to give you a temporary directory and cleanup when done.
    /// </summary>
    public class TemporaryDirectoryFixture : DirectoryLocation, IDisposable
    {
        public TemporaryDirectoryFixture()
            : base(Path.Combine(Path.GetTempPath(), Path.GetRandomFileName()))
        {
            Create();
        }

        public IFileLocation CreateFile(string filename)
        {
            var file = GetFileLocation(filename);
            file.Touch();
            return file;
        }

        public void Dispose()
        {
            Delete();
        }

        /// Override these setters to do nothing because AutoFixture will call them with garbage data
        public override FileAttributes Attributes
        {
            get { return base.Attributes; }
            set { }
        }

        public override DateTime CreationTime
        {
            get { return base.CreationTime; }
            set { }
        }

        public override DateTime CreationTimeUtc
        {
            get { return base.CreationTimeUtc; }
            set { }
        }

        public override DateTime LastAccessTime
        {
            get { return base.LastAccessTime; }
            set { }
        }

        public override DateTime LastAccessTimeUtc
        {
            get { return base.LastAccessTimeUtc; }
            set { }
        }

        public override DateTime LastWriteTime
        {
            get { return base.LastWriteTime; }
            set { }
        }

        public override DateTime LastWriteTimeUtc
        {
            get { return base.LastAccessTimeUtc; }
            set { }
        }
    }
}
