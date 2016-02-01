using System;
using System.IO;

namespace Isas.Shared
{
    /// <summary>
    /// Provides the base class for both FileLocation and DirectoryLocation objects
    /// </summary>
    public abstract class FileSystemLocationBase
    {
        /// <summary>
        /// Deletes a file or directory.
        /// </summary>
        public abstract void Delete();

        /// <summary>
        /// Refreshes the state of the object
        /// </summary>
        public abstract void Refresh();

        /// <summary>
        /// Gets or sets the attributes for the current file or directory.
        /// </summary>
        /// <value>
        /// FileAttributes of the current FileSystemInfo.
        /// </value>
        public abstract FileAttributes Attributes { get; set; }

        /// <summary>
        /// Gets or sets the creation time of the current file or directory.
        /// </summary>
        /// <value>
        /// The creation date and time of the current FileSystemInfo object.
        /// </value>
        public abstract DateTime CreationTime { get; set; }

        /// <summary>
        /// Gets or sets the creation time, in coordinated universal time (UTC), of the current file or directory.
        /// </summary>
        /// <value>
        /// The creation date and time in UTC format of the current FileSystemInfo object.
        /// </value>
        public abstract DateTime CreationTimeUtc { get; set; }

        /// <summary>
        /// Gets a value indicating whether the file or directory exists.
        /// </summary>
        /// <value>
        ///   <c>true</c> if the file or directory exists; otherwise, <c>false</c>.
        /// </value>
        public abstract bool Exists { get; }

        /// <summary>
        /// Gets the string representing the extension part of the file.
        /// </summary>
        /// <value>
        /// The Extension property returns the FileSystemInfo extension, including the period (.). For example, for a file c:\NewFile.txt, this property returns ".txt".
        /// </value>
        public abstract string Extension { get; }

        /// <summary>
        /// Gets the full path of the directory or file.
        /// </summary>
        /// <value>
        /// A string containing the full path.
        /// </value>
        public abstract string FullName { get; }

        /// <summary>
        /// Gets or sets the time the current file or directory was last accessed.
        /// </summary>
        /// <value>
        /// The time that the current file or directory was last accessed.
        /// </value>
        public abstract DateTime LastAccessTime { get; set; }

        /// <summary>
        /// Gets or sets the time, in coordinated universal time (UTC), that the current file or directory was last accessed.
        /// </summary>
        /// <value>
        /// The UTC time that the current file or directory was last accessed.
        /// </value>
        public abstract DateTime LastAccessTimeUtc { get; set; }

        /// <summary>
        /// Gets or sets the time when the current file or directory was last written to.
        /// </summary>
        /// <value>
        /// The time the current file was last written.
        /// </value>
        public abstract DateTime LastWriteTime { get; set; }

        /// <summary>
        /// Gets or sets the time, in coordinated universal time (UTC), when the current file or directory was last written to.
        /// </summary>
        /// <value>
        /// The UTC time when the current file was last written to.
        /// </value>
        public abstract DateTime LastWriteTimeUtc { get; set; }

        /// <summary>
        /// For files, gets the name of the file. For directories, gets the name of the last directory in the hierarchy if a hierarchy exists. Otherwise, the <c>Name</c> property gets the name of the directory.
        /// </summary>
        /// <value>
        /// For a directory, <c>Name</c> returns only the name of the parent directory, such as Dir, not c:\Dir. For a subdirectory, <c>Name</c> returns only the name of the subdirectory, such as Sub1, not c:\Dir\Sub1.
        /// For a file, <c>Name</c> returns only the file name and file name extension, such as MyFile.txt, not c:\Dir\Myfile.txt.
        /// </value>
        public abstract string Name { get; }
    }
}
