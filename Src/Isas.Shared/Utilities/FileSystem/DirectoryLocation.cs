using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Security.AccessControl;
using System.Threading;

namespace Isas.Shared
{
    /// <summary>
    /// Interface for exposing instance methods for creating, moving, and enumerating through directories and subdirectories.
    /// </summary>
    public interface IDirectoryLocation
    {
        /// <summary>
        /// Deletes this IDirectoryLocation if it is empty.
        /// </summary>
        void Delete();

        /// <summary>
        /// Gets or sets the attributes for the current file or directory.
        /// </summary>
        /// <value>
        /// FileAttributes of the current FileSystemInfo.
        /// </value>
        FileAttributes Attributes { get; set; }

        /// <summary>
        /// Gets or sets the creation time of the current file or directory.
        /// </summary>
        /// <value>
        /// The creation date and time of the current FileSystemLocationBase object.
        /// </value>
        DateTime CreationTime { get; set; }

        /// <summary>
        /// Gets or sets the creation time, in coordinated universal time (UTC), of the current file or directory.
        /// </summary>
        /// <value>
        /// The creation date and time in UTC format of the current FileSystemLocationBase object.
        /// </value>
        DateTime CreationTimeUtc { get; set; }

        /// <summary>
        /// Gets a value indicating whether the directory exists.
        /// </summary>
        /// <value>
        ///   <c>true</c> if the directory exists; otherwise, <c>false</c>.
        /// </value>
        bool Exists { get; }

        /// <summary>
        /// Gets the string representing the extension part of the file.
        /// </summary>
        /// <value>
        /// A string containing the FileSystemLocationBase extension.
        /// </value>
        string Extension { get; }

        /// <summary>
        /// Gets the full path of the directory or file.
        /// </summary>
        /// <value>
        /// A string containing the full path.
        /// </value>
        string FullName { get; }

        /// <summary>
        /// Gets or sets the time the current file or directory was last accessed.
        /// </summary>
        /// <value>
        /// The time that the current file or directory was last accessed.
        /// </value>
        DateTime LastAccessTime { get; set; }

        /// <summary>
        /// Gets or sets the time, in coordinated universal time (UTC), that the current file or directory was last accessed.
        /// </summary>
        /// <value>
        /// The UTC time that the current file or directory was last accessed.
        /// </value>
        DateTime LastAccessTimeUtc { get; set; }

        /// <summary>
        /// Gets or sets the time when the current file or directory was last written to.
        /// </summary>
        /// <value>
        /// The time the current file was last written.
        /// </value>
        DateTime LastWriteTime { get; set; }

        /// <summary>
        /// Gets or sets the time, in coordinated universal time (UTC), when the current file or directory was last written to.
        /// </summary>
        /// <value>
        /// The UTC time when the current file was last written to.
        /// </value>
        DateTime LastWriteTimeUtc { get; set; }

        /// <summary>
        /// Gets the name of this IDirectoryLocation instance.
        /// </summary>
        /// <value>
        /// The directory name.
        /// </value>
        string Name { get; }

        /// <summary>
        /// Gets the parent directory of a specified subdirectory.
        /// </summary>
        /// <value>
        /// The parent directory, or null if the path is <c>null</c> or if the file path denotes a root (such as "\", "C:", or * "\\server\share").
        /// </value>
        IDirectoryLocation Parent { get; }

        /// <summary>
        /// Gets the root portion of the directory.
        /// </summary>
        /// <value>
        /// An object that represents the root of the directory.
        /// </value>
        IDirectoryLocation Root { get; }

        /// <summary>
        /// Creates a directory. If the directory already exists, this method does nothing.
        /// </summary>
        void Create();

        /// <summary>
        /// Moves a directory into a parent directory. If paths are not on the same filesystem volume this is equivalent to CopyInto() then Delete()
        /// Any exisiting directory inside the parentDirectoryLocation with the same name is first removed
        /// Returns the directory in its new location
        /// </summary>
        IDirectoryLocation MoveInto(IDirectoryLocation parentDirectoryLocation);

        /// <summary>
        /// Moves a directory to a new location. If paths are not on the same filesystem volume this is equivalent to CopyTo() then Delete()
        /// Input is the directory in its new location and not the parent directory
        /// Any exisiting newDirectoryLocation is first removed
        /// Returns newDirectoryLocation to enable chaining
        /// </summary>
        IDirectoryLocation MoveTo(IDirectoryLocation newDirectoryLocation);

        /// <summary>
        /// Copy a directory into some parent directory. 
        /// Any exisiting directory inside the parentDirectoryLocation with the same name is first removed
        /// Returns the copied directory
        /// </summary>
        IDirectoryLocation CopyInto(IDirectoryLocation parentDirectoryLocation);

        /// <summary>
        /// Copy a directory to a new location. Input is the directory location for the resulting copied directory and not the parent directory. 
        /// Any exisiting newDirectoryLocation is first removed
        /// Returns newDirectoryLocation to enable chaining
        /// </summary>
        IDirectoryLocation CopyTo(IDirectoryLocation newDirectoryLocation);

        /// <summary>
        /// Creates a directory using a DirectorySecurity object.
        /// </summary>
        /// <param name="directorySecurity">The access control to apply to the directory.</param>
        void Create(DirectorySecurity directorySecurity);

        /// <summary>
        /// Creates a subdirectory or subdirectories on the specified path. The specified path can be relative to this instance of the IDirectoryLocation interface.
        /// </summary>
        /// <param name="name">The specified path. This cannot be a different disk volume or Universal Naming Convention (UNC) name.</param>
        /// <returns>The last directory specified in <paramref name="name"/>.</returns>
        IDirectoryLocation CreateSubdirectory(string name);

        /// <summary>
        /// Creates a subdirectory or subdirectories on the specified path with the specified security. The specified path can be relative to this instance of the IDirectoryLocation interface.
        /// </summary>
        /// <param name="path">The specified path. This cannot be a different disk volume or Universal Naming Convention (UNC) name.</param>
        /// <param name="directorySecurity">The security to apply.</param>
        /// <returns>The last directory specified in <paramref name="path"/>.</returns>
        IDirectoryLocation CreateSubdirectory(string path, DirectorySecurity directorySecurity);

        /// <summary>
        /// Returns an IFileLocation representing a file named <paramref name="filename"/> in the directory of this instance.
        /// The IFileLocation is instantiated, but the underlying file is not created. 
        /// </summary>
        /// <param name="filename">The name of the file to be created.</param>
        /// <returns>An IFileLocation representing the file created.</returns>
        IFileLocation GetFileLocation(string filename);

        /// <summary>
        /// Returns an IDirectoryLocation representing a file named <paramref name="dirname"/> in the directory of this instance.
        /// The IDirectoryLocation is instantiated, but the underlying directory is not created. 
        /// </summary>
        /// <param name="dirname">The name of the directory to be created.</param>
        /// <returns>An IDirectoryLocation representing the directory path.</returns>
        IDirectoryLocation GetDirectoryLocation(string dirname);

        /// <summary>
        /// Returns an enumerable collection of directory information in the current directory.
        /// </summary>
        /// <returns>An enumerable collection of directories in the current directory.</returns>
        IEnumerable<IDirectoryLocation> EnumerateDirectories();

        /// <summary>
        /// Returns an enumerable collection of directory information that matches a specified search pattern.
        /// </summary>
        /// <param name="searchPattern">The search string to match against the names of directories. This parameter can contain a combination of valid literal path and wildcard (* and ?) characters, but doesn't support regular expressions.</param>
        /// <returns>An enumerable collection of directories that matches <paramref name="searchPattern"/>.</returns>
        IEnumerable<IDirectoryLocation> EnumerateDirectories(string searchPattern);

        /// <summary>
        /// Returns an enumerable collection of directory information that matches a specified search pattern and search subdirectory option.
        /// </summary>
        /// <param name="searchPattern">The search string to match against the names of directories. This parameter can contain a combination of valid literal path and wildcard (* and ?) characters, but doesn't support regular expressions.</param>
        /// <param name="searchOption">One of the enumeration values that specifies whether the search operation should include only the current directory or all subdirectories.</param>
        /// <returns></returns>
        IEnumerable<IDirectoryLocation> EnumerateDirectories(string searchPattern, SearchOption searchOption);

        /// <summary>
        /// Returns an enumerable collection of file information in the current directory.
        /// </summary>
        /// <returns>An enumerable collection of the files in the current directory.</returns>
        IEnumerable<IFileLocation> EnumerateFiles();

        /// <summary>
        /// Returns an enumerable collection of file information that matches a search pattern.
        /// </summary>
        /// <param name="searchPattern">The search string to match against the names of files. This parameter can contain a combination of valid literal path and wildcard (* and ?) characters, but doesn't support regular expressions.</param>
        /// <returns>An enumerable collection of files that matches <paramref name="searchPattern"/>.</returns>
        IEnumerable<IFileLocation> EnumerateFiles(string searchPattern);

        /// <summary>
        /// Returns an enumerable collection of file information that matches a specified search pattern and search subdirectory option.
        /// </summary>
        /// <param name="searchPattern">The search string to match against the names of files. This parameter can contain a combination of valid literal path and wildcard (* and ?) characters, but doesn't support regular expressions.</param>
        /// <param name="searchOption">One of the enumeration values that specifies whether the search operation should include only the current directory or all subdirectories. </param>
        /// <returns>An enumerable collection of files that matches searchPattern and searchOption.</returns>
        IEnumerable<IFileLocation> EnumerateFiles(string searchPattern, SearchOption searchOption);

        /// <summary>
        /// Returns an enumerable collection of file system information in the current directory.
        /// </summary>
        /// <returns>An enumerable collection of file system information in the current directory. </returns>
        IEnumerable<FileSystemLocationBase> EnumerateFileSystemInfos();

        /// <summary>
        /// Returns an enumerable collection of file system information that matches a specified search pattern.
        /// </summary>
        /// <param name="searchPattern">The search string to match against the names of directories. This parameter can contain a combination of valid literal path and wildcard (* and ?) characters (see Remarks), but doesn't support regular expressions.</param>
        /// <returns>An enumerable collection of file system information objects that matches <paramref name="searchPattern"/>.</returns>
        IEnumerable<FileSystemLocationBase> EnumerateFileSystemInfos(string searchPattern);

        /// <summary>
        /// Returns an enumerable collection of file system information that matches a specified search pattern and search subdirectory option.
        /// </summary>
        /// <param name="searchPattern">The search string to match against the names of directories. This parameter can contain a combination of valid literal path and wildcard (* and ?) characters (see Remarks), but doesn't support regular expressions.</param>
        /// <param name="searchOption">One of the enumeration values that specifies whether the search operation should include only the current directory or all subdirectories.</param>
        /// <returns>An enumerable collection of file system information objects that matches <paramref name="searchPattern"/> and <paramref name="searchOption"/>.</returns>
        IEnumerable<FileSystemLocationBase> EnumerateFileSystemInfos(string searchPattern, SearchOption searchOption);

        /// <summary>
        /// Gets a DirectorySecurity object that encapsulates the access control list (ACL) entries for the directory described by the current DirectoryInfo object.
        /// </summary>
        /// <returns>A DirectorySecurity object that encapsulates the access control rules for the directory.</returns>
        DirectorySecurity GetAccessControl();

        /// <summary>
        /// Gets a DirectorySecurity object that encapsulates the specified type of access control list (ACL) entries for the directory described by the current DirectoryInfo object.
        /// </summary>
        /// <param name="includeSections">One of the AccessControlSections values that specifies the type of access control list (ACL) information to receive.</param>
        /// <returns>A DirectorySecurity object that encapsulates the access control rules for the file.</returns>
        DirectorySecurity GetAccessControl(AccessControlSections includeSections);

        /// <summary>
        /// Returns the subdirectories of the current directory.
        /// </summary>
        /// <returns>An array of IDirectoryLocation objects.</returns>
        IDirectoryLocation[] GetDirectories();

        /// <summary>
        /// Returns an array of directories in the current DirectoryInfo matching the given search criteria.
        /// </summary>
        /// <param name="searchPattern">The search string to match against the names of directories. This parameter can contain a combination of valid literal path and wildcard (* and ?) characters, but doesn't support regular expressions.</param>
        /// <returns>An array of type <c>IDirectoryLocation</c> matching searchPattern.</returns>
        IDirectoryLocation[] GetDirectories(string searchPattern);

        /// <summary>
        /// Returns an array of directories in the current IDirectoryLocation matching the given search criteria and using a value to determine whether to search subdirectories.
        /// </summary>
        /// <param name="searchPattern">The search string to match against the names of directories. This parameter can contain a combination of valid literal path and wildcard (* and ?) characters, but doesn't support regular expressions.</param>
        /// <param name="searchOption">One of the enumeration values that specifies whether the search operation should include only the current directory or all subdirectories.</param>
        /// <returns>An array of type <c>IDirectoryLocation</c> matching searchPattern.</returns>
        IDirectoryLocation[] GetDirectories(string searchPattern, SearchOption searchOption);

        /// <summary>
        /// Returns a file list from the current directory.
        /// </summary>
        /// <returns>An array of type IFileLocation.</returns>
        IFileLocation[] GetFiles();

        /// <summary>
        /// Returns a file list from the current directory matching the given search pattern.
        /// </summary>
        /// <param name="searchPattern">The search string to match against the names of files. This parameter can contain a combination of valid literal path and wildcard (* and ?) characters, but doesn't support regular expressions.</param>
        /// <returns>An array of type IFileLocation.</returns>
        IFileLocation[] GetFiles(string searchPattern);

        /// <summary>
        /// Returns a file list from the current directory matching the given search pattern and using a value to determine whether to search subdirectories.
        /// </summary>
        /// <param name="searchPattern">The search string to match against the names of files. This parameter can contain a combination of valid literal path and wildcard (* and ?) characters, but doesn't support regular expressions.</param>
        /// <param name="searchOption">One of the enumeration values that specifies whether the search operation should include only the current directory or all subdirectories.</param>
        /// <returns>An array of type IFileLocation.</returns>
        IFileLocation[] GetFiles(string searchPattern, SearchOption searchOption);

        /// <summary>
        /// Returns an array of strongly typed FileSystemLocationBase entries representing all the files and subdirectories in a directory.
        /// </summary>
        /// <returns>An array of strongly typed FileSystemLocationBase entries.</returns>
        FileSystemLocationBase[] GetFileSystemInfos();

        /// <summary>
        /// Retrieves an array of strongly typed FileSystemLocationBase objects representing the files and subdirectories that match the specified search criteria.
        /// </summary>
        /// <param name="searchPattern">The search string to match against the names of directories and files. This parameter can contain a combination of valid literal path and wildcard (* and ?) characters, but doesn't support regular expressions.</param>
        /// <returns>An array of strongly typed FileSystemLocationBase objects matching the search criteria.</returns>
        FileSystemLocationBase[] GetFileSystemInfos(string searchPattern);

        /// <summary>
        /// Retrieves an array of FileSystemLocationBase objects that represent the files and subdirectories matching the specified search criteria.
        /// </summary>
        /// <param name="searchPattern">The search string to match against the names of directories and files. This parameter can contain a combination of valid literal path and wildcard (* and ?) characters, but doesn't support regular expressions.</param>
        /// <param name="searchOption">One of the enumeration values that specifies whether the search operation should include only the current directory or all subdirectories. The default value is TopDirectoryOnly.</param>
        /// <returns>An array of file system entries that match the search criteria.</returns>
        FileSystemLocationBase[] GetFileSystemInfos(string searchPattern, SearchOption searchOption);

        /// <summary>
        /// Applies access control list (ACL) entries described by a DirectorySecurity object to the directory described by the current IDirectoryLocation object.
        /// </summary>
        /// <param name="directorySecurity">An object that describes an ACL entry to apply to the directory.</param>
        void SetAccessControl(DirectorySecurity directorySecurity);

        /* Requires the Mono.Unix assembly
		/// <summary>
		/// Symlinks an IDirectoryLocation instance and its contents to a new path.
		/// </summary>
		/// <param name="destDirName">The name and path to which to symlink this directory.</param>
		/// <returns>A symlinked directory</returns>
		IDirectoryLocation SymlinkTo(string destDirName);
		*/

        /// <summary>
        /// Symlinks an existing file to a new location.
        /// If the symlink cannot be created (Windows may require elevated privileges to do this),
        ///  then a copy is created instead if <paramref name="copyOnFail"/> is true.
        /// </summary>
        /// <param name="sourceFile">The instance of IDirectoryLocation this method extends</param>
        /// <param name="destinationFolder">The IDirectoryLocation representing the symlink, 
        /// which points to <paramref name="sourceFile"/></param>
        /// <param name="copyOnFail">If true, copies the file if a symlink cannot be created.</param>
        /// <returns>An IDirectoryLocation of the destination file, or null if the symlink failed and <paramref name="copyOnFail"/> is false.</returns>
        IDirectoryLocation SymlinkTo(IDirectoryLocation destinationFolder, bool copyOnFail = true);

    }

    public static class DirectoryLocationExtensions
    {
        /// <summary>
        /// Returns a file location by combining all <paramref name="args"/> into a path relative
        ///  to this instance of IDirectoryLocation.
        /// </summary>
        /// <param name="dir"></param>
        /// <param name="args">The tokens of the relative file path. The last token represents the file name.</param>
        /// <returns></returns>
        public static IFileLocation GetFileLocation(this IDirectoryLocation dir, params string[] tokens)
        {
            IDirectoryLocation tempDir = dir;
            for (int pos = 0; pos < tokens.Count() - 1; pos++)
            {
                dir = dir.GetDirectoryLocation(tokens[pos]);
            }
            return dir.GetFileLocation(tokens.Last());
        }

        /// <summary>
        /// Returns a directory location by combining all <paramref name="args"/> into a path relative
        ///  to this instance of IDirectoryLocation. Returns only the path - the directory is not
        ///  guaranteed to exist.
        /// </summary>
        /// <param name="dir"></param>
        /// <param name="args">The tokens of the relative file path. The last token represents the file name.</param>
        /// <returns></returns>
        public static IDirectoryLocation GetDirectoryLocation(this IDirectoryLocation dir, params string[] tokens)
        {
            IDirectoryLocation tempDir = dir;
            foreach (string token in tokens)
            {
                tempDir = tempDir.GetDirectoryLocation(token);
            }
            return tempDir;
        }

        /// <summary>
        /// Returns a directory location by combining all <paramref name="args"/> into a path relative
        ///  to this instance of IDirectoryLocation. Creates the directory.
        /// </summary>
        /// <param name="dir"></param>
        /// <param name="tokens">The tokens of the relative file path. The last token represents the file name.</param>
        /// <returns></returns>
        public static IDirectoryLocation CreateSubdirectory(this IDirectoryLocation dir, params string[] tokens)
        {
            IDirectoryLocation tempDir = dir;
            foreach (string token in tokens)
            {
                tempDir = tempDir.CreateSubdirectory(token);
            }
            return tempDir;
        }

        /// <summary>
        /// Returns a directory location by combining all <paramref name="args"/> into a path relative
        ///  to this instance of IDirectoryLocation. Deletes any existing copy, then creates the directory.
        /// </summary>
        /// <param name="dir"></param>
        /// <param name="tokens">The tokens of the relative file path. The last token represents the file name.</param>
        /// <returns></returns>
        public static IDirectoryLocation CreateEmptySubdirectory(this IDirectoryLocation dir, params string[] tokens)
        {
            IDirectoryLocation tempDir = dir;
            foreach (string token in tokens)
            {
                tempDir = tempDir.GetDirectoryLocation(token).CreateClean();
            }
            return tempDir;
        }

        /// <summary>
        /// Creates a directory by first deleting any existing directory
        /// </summary>
        public static IDirectoryLocation CreateClean(this IDirectoryLocation dir)
        {
            dir.Delete();
            dir.Create();
            return dir;
        }
    }

    /// <summary>
    /// Exposes instance methods for creating, moving, and enumerating through directories and subdirectories.
    /// </summary>
    public class DirectoryLocation : FileSystemLocationBase, IDirectoryLocation
    {
        private readonly DirectoryInfo _instance;

        /// <summary>
        /// Initializes a new instance of the <see cref="DirectoryLocation"/> class.
        /// </summary>
        /// <param name="directoryInfo">A DirectoryInfo object to wrap.</param>
        public DirectoryLocation(DirectoryInfo directoryInfo)
        {
            _instance = directoryInfo;
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="DirectoryLocation"/> class.
        /// </summary>
        /// <param name="directoryPath">The directory path.</param>
        public DirectoryLocation(string directoryPath)
        {
            _instance = new DirectoryInfo(directoryPath);
        }

        /// <summary>
        /// Refreshes the state of the object.
        /// </summary>
        public override void Refresh()
        {
            _instance.Refresh();
        }

        /// <summary>
        /// Gets or sets the attributes for the current file or directory.
        /// </summary>
        /// <value>
        /// FileAttributes of the current FileSystemInfo.
        /// </value>
        public override FileAttributes Attributes
        {
            get { Refresh(); return _instance.Attributes; }
            set { _instance.Attributes = value; }
        }

        /// <summary>
        /// Gets or sets the creation time of the current file or directory.
        /// </summary>
        /// <value>
        /// The creation date and time of the current FileSystemLocationBase object.
        /// </value>
        public override DateTime CreationTime
        {
            get { Refresh(); return _instance.CreationTime; }
            set { _instance.CreationTime = value; }
        }

        /// <summary>
        /// Gets or sets the creation time, in coordinated universal time (UTC), of the current file or directory.
        /// </summary>
        /// <value>
        /// The creation date and time in UTC format of the current FileSystemLocationBase object.
        /// </value>
        public override DateTime CreationTimeUtc
        {
            get { Refresh(); return _instance.CreationTimeUtc; }
            set { _instance.CreationTimeUtc = value; }
        }

        /// <summary>
        /// Gets a value indicating whether the directory exists.
        /// </summary>
        /// <value>
        ///   <c>true</c> if the directory exists; otherwise, <c>false</c>.
        /// </value>
        public override bool Exists
        {
            get { Refresh(); return _instance.Exists; }
        }

        /// <summary>
        /// Gets the string representing the extension part of the file.
        /// </summary>
        /// <value>
        /// A string containing the FileSystemLocationBase extension.
        /// </value>
        public override string Extension
        {
            get { return _instance.Extension; }
        }

        /// <summary>
        /// Gets the full path of the directory or file.
        /// </summary>
        /// <value>
        /// A string containing the full path.
        /// </value>
        public override string FullName
        {
            get { return _instance.FullName; }
        }

        /// <summary>
        /// Gets or sets the time the current file or directory was last accessed.
        /// </summary>
        /// <value>
        /// The time that the current file or directory was last accessed.
        /// </value>
        public override DateTime LastAccessTime
        {
            get { Refresh(); return _instance.LastAccessTime; }
            set { _instance.LastAccessTime = value; }
        }

        /// <summary>
        /// Gets or sets the time, in coordinated universal time (UTC), that the current file or directory was last accessed.
        /// </summary>
        /// <value>
        /// The UTC time that the current file or directory was last accessed.
        /// </value>
        public override DateTime LastAccessTimeUtc
        {
            get { Refresh(); return _instance.LastAccessTimeUtc; }
            set { _instance.LastAccessTimeUtc = value; }
        }

        /// <summary>
        /// Gets or sets the time when the current file or directory was last written to.
        /// </summary>
        /// <value>
        /// The time the current file was last written.
        /// </value>
        public override DateTime LastWriteTime
        {
            get { Refresh(); return _instance.LastWriteTime; }
            set { _instance.LastWriteTime = value; }
        }

        /// <summary>
        /// Gets or sets the time, in coordinated universal time (UTC), when the current file or directory was last written to.
        /// </summary>
        /// <value>
        /// The UTC time when the current file was last written to.
        /// </value>
        public override DateTime LastWriteTimeUtc
        {
            get { Refresh(); return _instance.LastWriteTimeUtc; }
            set { _instance.LastWriteTimeUtc = value; }
        }

        /// <summary>
        /// Gets the name of this DirectoryLocation instance.
        /// </summary>
        /// <value>
        /// The directory name.
        /// </value>
        public override string Name
        {
            get { return _instance.Name; }
        }

        /// <summary>
        /// Deletes this DirectoryLocation including all files and subdirectories.
        /// </summary>
        public override void Delete()
        {
            if (Exists)
            {
                DeleteFilesAndFoldersBruteForce(FullName);
                //DeleteFilesAndFoldersRecursively(FullName);
                // NOTE: below does not always work on windows for some reason (see unit test: TestFileLocationMultipleTimes)
                //_instance.Delete(true);
            }
        }

        //// this method doesn't seem to want to delete symlinks on linux :(
        //private static void DeleteFilesAndFoldersRecursively(string target_dir)
        //{
        //    foreach (string file in Directory.GetFiles(target_dir))
        //    {
        //        File.Delete(file);
        //    }

        //    foreach (string subDir in Directory.GetDirectories(target_dir))
        //    {
        //        DeleteFilesAndFoldersRecursively(subDir);
        //    }

        //    Thread.Sleep(10); // This makes the difference between whether it works or not. Sleep(0) is not enough.
        //    Directory.Delete(target_dir);
        //}

        //this also seems to be faster than the recursive method above
        private static void DeleteFilesAndFoldersBruteForce(string targetDir)
        {
            int maxAttempt = 3;
            int attempt = 1;
            while (true)
            {
                try
                {
                    Directory.Delete(targetDir, true);
                    return;
                }
                catch (Exception ex)
                {
                    if (attempt == maxAttempt)
                        throw new ApplicationException($"Unable to recursively remove {targetDir} after {attempt} attempts. Do you have a circular symbolic link? Detailed error: {ex}", ex);
                    attempt++;
                    Thread.Sleep(10 * attempt * attempt);
                }
            }
        }

        /// <summary>
        /// Creates a directory. If the directory already exists, this method does nothing.
        /// </summary>
        public void Create()
        {
            _instance.Create();
        }

        /// <summary>
        /// Moves a directory into a parent directory. If paths are not on the same filesystem volume this is equivalent to CopyInto() then Delete()
        /// If the parentDirectoryLocation already contains a directory with the same it is first removed
        /// Returns the directory in its new location
        /// </summary>
        public IDirectoryLocation MoveInto(IDirectoryLocation parentDirectoryLocation)
        {
            return MoveTo(parentDirectoryLocation.CreateSubdirectory(Name));
        }

        /// <summary>
        /// Moves a directory to a new location. If paths are not on the same filesystem volume this is equivalent to CopyTo() then Delete()
        /// Input is the directory in its new location and not the parent directory
        /// If the newDirectoryLocation already exists it is first removed
        /// Returns newDirectoryLocation to enable chaining
        /// </summary>
        public IDirectoryLocation MoveTo(IDirectoryLocation newDirectoryLocation)
        {
            if (newDirectoryLocation.FullName == FullName) return this;
            if (newDirectoryLocation.Parent != null)
                newDirectoryLocation.Parent.Create();
            newDirectoryLocation.Delete();
            try
            {
                Directory.Move(FullName, newDirectoryLocation.FullName);
            }
            catch (IOException)
            {
                CopyTo(newDirectoryLocation);
                Delete();
            }
            return newDirectoryLocation;
        }

        /// <summary>
        /// Copy a directory into some parent directory. 
        /// If a directory with the same name exists in the parent directory it is first removed. 
        /// Returns the copied directory
        /// </summary>
        public IDirectoryLocation CopyInto(IDirectoryLocation parentDirectoryLocation)
        {
            return CopyTo(parentDirectoryLocation.CreateSubdirectory(Name));
        }

        /// <summary>
        /// Copy a directory to a new location. Input is the directory location for the resulting copied directory and not the parent directory. 
        /// If newDirectoryLocation already exists it is first removed. 
        /// Returns newDirectoryLocation to enable chaining
        /// </summary>
        public IDirectoryLocation CopyTo(IDirectoryLocation newDirectoryLocation)
        {
            if (newDirectoryLocation.FullName == FullName) return this;
            newDirectoryLocation.CreateClean();

            foreach (IFileLocation file in EnumerateFiles())
            {
                file.CopyTo(newDirectoryLocation.GetFileLocation(file.Name));
            }
            foreach (IDirectoryLocation subdir in EnumerateDirectories())
            {
                subdir.CopyInto(newDirectoryLocation);
            }
            return newDirectoryLocation;
        }

        /// <summary>
        /// Creates a directory using a DirectorySecurity object.
        /// </summary>
        /// <param name="directorySecurity">The access control to apply to the directory.</param>
        public void Create(DirectorySecurity directorySecurity)
        {
            _instance.Create(directorySecurity);
            Refresh();
        }

        /// <summary>
        /// Creates a subdirectory or subdirectories on the specified path. The specified path can be relative to this instance of the DirectoryLocation interface.
        /// </summary>
        /// <param name="path">The specified path. This cannot be a different disk volume or Universal Naming Convention (UNC) name.</param>
        /// <returns>
        /// The last directory specified in <paramref name="path" />.
        /// </returns>
        public IDirectoryLocation CreateSubdirectory(string path)
        {
            return new DirectoryLocation(_instance.CreateSubdirectory(path));
        }

        /// <summary>
        /// Creates a subdirectory or subdirectories on the specified path with the specified security. The specified path can be relative to this instance of the DirectoryLocation interface.
        /// </summary>
        /// <param name="path">The specified path. This cannot be a different disk volume or Universal Naming Convention (UNC) name.</param>
        /// <param name="directorySecurity">The security to apply.</param>
        /// <returns>
        /// The last directory specified in <paramref name="path" />.
        /// </returns>
        public IDirectoryLocation CreateSubdirectory(string path, DirectorySecurity directorySecurity)
        {
            return new DirectoryLocation(_instance.CreateSubdirectory(path, directorySecurity));
        }

        /// <summary>
        /// Returns an enumerable collection of directory information in the current directory.
        /// </summary>
        /// <returns>
        /// An enumerable collection of directories in the current directory.
        /// </returns>
        public IEnumerable<IDirectoryLocation> EnumerateDirectories()
        {
            return _instance.EnumerateDirectories().Select(directoryInfo => new DirectoryLocation(directoryInfo));
        }

        /// <summary>
        /// Returns an enumerable collection of directory information that matches a specified search pattern.
        /// </summary>
        /// <param name="searchPattern">The search string to match against the names of directories. This parameter can contain a combination of valid literal path and wildcard (* and ?) characters, but doesn't support regular expressions.</param>
        /// <returns>
        /// An enumerable collection of directories that matches <paramref name="searchPattern" />.
        /// </returns>
        public IEnumerable<IDirectoryLocation> EnumerateDirectories(string searchPattern)
        {
            return _instance.EnumerateDirectories(searchPattern).Select(directoryInfo => new DirectoryLocation(directoryInfo));
        }

        /// <summary>
        /// Returns an enumerable collection of directory information that matches a specified search pattern and search subdirectory option.
        /// </summary>
        /// <param name="searchPattern">The search string to match against the names of directories. This parameter can contain a combination of valid literal path and wildcard (* and ?) characters, but doesn't support regular expressions.</param>
        /// <param name="searchOption">One of the enumeration values that specifies whether the search operation should include only the current directory or all subdirectories.</param>
        /// <returns></returns>
        public IEnumerable<IDirectoryLocation> EnumerateDirectories(string searchPattern, SearchOption searchOption)
        {
            return _instance.EnumerateDirectories(searchPattern, searchOption).Select(directoryInfo => new DirectoryLocation(directoryInfo));
        }

        /// <summary>
        /// Returns an enumerable collection of file information in the current directory.
        /// </summary>
        /// <returns>
        /// An enumerable collection of the files in the current directory.
        /// </returns>
        public IEnumerable<IFileLocation> EnumerateFiles()
        {
            return _instance.EnumerateFiles().Select(fileInfo => new FileLocation(fileInfo));
        }

        /// <summary>
        /// Returns an enumerable collection of file information that matches a search pattern.
        /// </summary>
        /// <param name="searchPattern">The search string to match against the names of files. This parameter can contain a combination of valid literal path and wildcard (* and ?) characters, but doesn't support regular expressions.</param>
        /// <returns>
        /// An enumerable collection of files that matches <paramref name="searchPattern" />.
        /// </returns>
        public IEnumerable<IFileLocation> EnumerateFiles(string searchPattern)
        {
            return _instance.EnumerateFiles(searchPattern).Select(fileInfo => new FileLocation(fileInfo));
        }

        /// <summary>
        /// Returns an enumerable collection of file information that matches a specified search pattern and search subdirectory option.
        /// </summary>
        /// <param name="searchPattern">The search string to match against the names of files. This parameter can contain a combination of valid literal path and wildcard (* and ?) characters, but doesn't support regular expressions.</param>
        /// <param name="searchOption">One of the enumeration values that specifies whether the search operation should include only the current directory or all subdirectories.</param>
        /// <returns>
        /// An enumerable collection of files that matches searchPattern and searchOption.
        /// </returns>
        public IEnumerable<IFileLocation> EnumerateFiles(string searchPattern, SearchOption searchOption)
        {
            return _instance.EnumerateFiles(searchPattern, searchOption).Select(fileInfo => new FileLocation(fileInfo));
        }

        /// <summary>
        /// Returns an enumerable collection of file system information in the current directory.
        /// </summary>
        /// <returns>
        /// An enumerable collection of file system information in the current directory.
        /// </returns>
        public IEnumerable<FileSystemLocationBase> EnumerateFileSystemInfos()
        {
            return _instance.EnumerateFileSystemInfos().WrapFileSystemInfos();
        }

        /// <summary>
        /// Returns an enumerable collection of file system information that matches a specified search pattern.
        /// </summary>
        /// <param name="searchPattern">The search string to match against the names of directories. This parameter can contain a combination of valid literal path and wildcard (* and ?) characters (see Remarks), but doesn't support regular expressions.</param>
        /// <returns>
        /// An enumerable collection of file system information objects that matches <paramref name="searchPattern" />.
        /// </returns>
        public IEnumerable<FileSystemLocationBase> EnumerateFileSystemInfos(string searchPattern)
        {
            return _instance.EnumerateFileSystemInfos(searchPattern).WrapFileSystemInfos();
        }

        /// <summary>
        /// Returns an enumerable collection of file system information that matches a specified search pattern and search subdirectory option.
        /// </summary>
        /// <param name="searchPattern">The search string to match against the names of directories. This parameter can contain a combination of valid literal path and wildcard (* and ?) characters (see Remarks), but doesn't support regular expressions.</param>
        /// <param name="searchOption">One of the enumeration values that specifies whether the search operation should include only the current directory or all subdirectories.</param>
        /// <returns>
        /// An enumerable collection of file system information objects that matches <paramref name="searchPattern" /> and <paramref name="searchOption" />.
        /// </returns>
        public IEnumerable<FileSystemLocationBase> EnumerateFileSystemInfos(string searchPattern, SearchOption searchOption)
        {
            return _instance.EnumerateFileSystemInfos(searchPattern, searchOption).WrapFileSystemInfos();
        }

        /// <summary>
        /// Gets a DirectorySecurity object that encapsulates the access control list (ACL) entries for the directory described by the current DirectoryInfo object.
        /// </summary>
        /// <returns>
        /// A DirectorySecurity object that encapsulates the access control rules for the directory.
        /// </returns>
        public DirectorySecurity GetAccessControl()
        {
            return _instance.GetAccessControl();
        }

        /// <summary>
        /// Gets a DirectorySecurity object that encapsulates the specified type of access control list (ACL) entries for the directory described by the current DirectoryInfo object.
        /// </summary>
        /// <param name="includeSections">One of the AccessControlSections values that specifies the type of access control list (ACL) information to receive.</param>
        /// <returns>
        /// A DirectorySecurity object that encapsulates the access control rules for the file.
        /// </returns>
        public DirectorySecurity GetAccessControl(AccessControlSections includeSections)
        {
            return _instance.GetAccessControl(includeSections);
        }

        /// <summary>
        /// Returns the subdirectories of the current directory.
        /// </summary>
        /// <returns>
        /// An array of IDirectoryLocation objects.
        /// </returns>
        public IDirectoryLocation[] GetDirectories()
        {
            return _instance.GetDirectories().WrapDirectories();
        }

        /// <summary>
        /// Returns an array of directories in the current DirectoryInfo matching the given search criteria.
        /// </summary>
        /// <param name="searchPattern">The search string to match against the names of directories. This parameter can contain a combination of valid literal path and wildcard (* and ?) characters, but doesn't support regular expressions.</param>
        /// <returns>
        /// An array of type <c>IDirectoryLocation</c> matching searchPattern.
        /// </returns>
        public IDirectoryLocation[] GetDirectories(string searchPattern)
        {
            return _instance.GetDirectories(searchPattern).WrapDirectories();
        }

        /// <summary>
        /// Returns an array of directories in the current DirectoryLocation matching the given search criteria and using a value to determine whether to search subdirectories.
        /// </summary>
        /// <param name="searchPattern">The search string to match against the names of directories. This parameter can contain a combination of valid literal path and wildcard (* and ?) characters, but doesn't support regular expressions.</param>
        /// <param name="searchOption">One of the enumeration values that specifies whether the search operation should include only the current directory or all subdirectories.</param>
        /// <returns>
        /// An array of type <c>IDirectoryLocation</c> matching searchPattern.
        /// </returns>
        public IDirectoryLocation[] GetDirectories(string searchPattern, SearchOption searchOption)
        {
            return _instance.GetDirectories(searchPattern, searchOption).WrapDirectories();
        }

        /// <summary>
        /// Returns a file list from the current directory.
        /// </summary>
        /// <returns>
        /// An array of type IFileLocation.
        /// </returns>
        public IFileLocation[] GetFiles()
        {
            return _instance.GetFiles().WrapFiles();
        }

        /// <summary>
        /// Returns a file list from the current directory matching the given search pattern.
        /// </summary>
        /// <param name="searchPattern">The search string to match against the names of files. This parameter can contain a combination of valid literal path and wildcard (* and ?) characters, but doesn't support regular expressions.</param>
        /// <returns>
        /// An array of type IFileLocation.
        /// </returns>
        public IFileLocation[] GetFiles(string searchPattern)
        {
            return _instance.GetFiles(searchPattern).WrapFiles();
        }

        /// <summary>
        /// Returns a file list from the current directory matching the given search pattern and using a value to determine whether to search subdirectories.
        /// </summary>
        /// <param name="searchPattern">The search string to match against the names of files. This parameter can contain a combination of valid literal path and wildcard (* and ?) characters, but doesn't support regular expressions.</param>
        /// <param name="searchOption">One of the enumeration values that specifies whether the search operation should include only the current directory or all subdirectories.</param>
        /// <returns>
        /// An array of type IFileLocation.
        /// </returns>
        public IFileLocation[] GetFiles(string searchPattern, SearchOption searchOption)
        {
            return _instance.GetFiles(searchPattern, searchOption).WrapFiles();
        }

        /// <summary>
        /// Returns an array of strongly typed FileSystemLocationBase entries representing all the files and subdirectories in a directory.
        /// </summary>
        /// <returns>
        /// An array of strongly typed FileSystemLocationBase entries.
        /// </returns>
        public FileSystemLocationBase[] GetFileSystemInfos()
        {
            return _instance.GetFileSystemInfos().WrapFileSystemInfos();
        }

        /// <summary>
        /// Retrieves an array of strongly typed FileSystemLocationBase objects representing the files and subdirectories that match the specified search criteria.
        /// </summary>
        /// <param name="searchPattern">The search string to match against the names of directories and files. This parameter can contain a combination of valid literal path and wildcard (* and ?) characters, but doesn't support regular expressions.</param>
        /// <returns>
        /// An array of strongly typed FileSystemLocationBase objects matching the search criteria.
        /// </returns>
        public FileSystemLocationBase[] GetFileSystemInfos(string searchPattern)
        {
            return _instance.GetFileSystemInfos(searchPattern).WrapFileSystemInfos();
        }

        /// <summary>
        /// Retrieves an array of FileSystemLocationBase objects that represent the files and subdirectories matching the specified search criteria.
        /// </summary>
        /// <param name="searchPattern">The search string to match against the names of directories and files. This parameter can contain a combination of valid literal path and wildcard (* and ?) characters, but doesn't support regular expressions.</param>
        /// <param name="searchOption">One of the enumeration values that specifies whether the search operation should include only the current directory or all subdirectories. The default value is TopDirectoryOnly.</param>
        /// <returns>
        /// An array of file system entries that match the search criteria.
        /// </returns>
        public FileSystemLocationBase[] GetFileSystemInfos(string searchPattern, SearchOption searchOption)
        {
            return _instance.GetFileSystemInfos(searchPattern, searchOption).WrapFileSystemInfos();
        }

        /// <summary>
        /// Applies access control list (ACL) entries described by a DirectorySecurity object to the directory described by the current DirectoryLocation object.
        /// </summary>
        /// <param name="directorySecurity">An object that describes an ACL entry to apply to the directory</param>
        public void SetAccessControl(DirectorySecurity directorySecurity)
        {
            _instance.SetAccessControl(directorySecurity);
        }

        /// <summary>
        /// Gets the parent directory of a specified subdirectory.
        /// </summary>
        /// <value>
        /// The parent directory, or null if the path is <c>null</c> or if the file path denotes a root (such as "\", "C:", or * "\\server\share").
        /// </value>
        public IDirectoryLocation Parent
        {
            get { return new DirectoryLocation(_instance.Parent); }
        }

        /// <summary>
        /// Gets the root portion of the directory.
        /// </summary>
        /// <value>
        /// An object that represents the root of the directory.
        /// </value>
        public IDirectoryLocation Root
        {
            get { return new DirectoryLocation(_instance.Root); }
        }

        /* Requires Mono.Unix assembly
		[DllImport("kernel32.dll")]
		private static extern bool CreateSymbolicLink(string lpSymlinkFileName, string lpTargetFileName, int dwFlags);
		/// <summary>
		/// Symlinks an DirectoryLocation instance and its contents to a new path.
		/// </summary>
		/// <param name="destDirName">The name and path to which to symlink this directory.</param>
		/// <returns>
		/// A symlinked directory
		/// </returns>
		public IDirectoryLocation SymlinkTo(string destDirName)
		{
			if (Type.GetType("Mono.Runtime") != null)
			{
				var directory = new UnixDirectoryInfo(_instance.FullName).CreateSymbolicLink(destDirName);
				return new DirectoryLocation(directory.FullName);
			}

			CreateSymbolicLink(destDirName, _instance.FullName, 1);
			return new DirectoryLocation(destDirName);
		}
		*/

        /// <summary>
        /// Symlinks an existing file to a new location.
        /// If the symlink cannot be created (Windows may require elevated privileges to do this),
        ///  then a copy is created instead if <paramref name="copyOnFail"/> is true.
        /// </summary>
        /// <param name="sourceFile">The instance of IDirectoryLocation this method extends</param>
        /// <param name="destinationFolder">The IDirectoryLocation representing the symlink, 
        /// which points to <paramref name="sourceFile"/></param>
        /// <param name="copyOnFail">If true, copies the file if a symlink cannot be created.</param>
        /// <returns>An IDirectoryLocation of the destination file, or null if the symlink failed and <paramref name="copyOnFail"/> is false.</returns>
        public IDirectoryLocation SymlinkTo(IDirectoryLocation destinationFolder, bool copyOnFail = true)
        {
            return (Illumina.SecondaryAnalysis.Utilities.CreateSymbolicLink(destinationFolder.FullName, FullName, copyOnFail))
                    ? destinationFolder : null;
        }

        /// <summary>
        /// Return the path to a file within the current instance.
        /// </summary>
        /// <param name="filename">The name of the file</param>
        public IFileLocation GetFileLocation(string filename)
        {
            return new FileLocation(Path.Combine(FullName, filename));
        }

        /// <summary>
        /// Return the path to a directory within the current instance.
        /// </summary>
        /// <param name="dirname">The name of the subdirectory</param>
        public IDirectoryLocation GetDirectoryLocation(string dirname)
        {
            return new DirectoryLocation(Path.Combine(FullName, dirname));
        }

        public override string ToString()
        {
            return FullName;
        }

        public override bool Equals(object obj)
        {
            if (typeof(IDirectoryLocation).IsAssignableFrom(obj.GetType()))
            {
                return FullName.Equals(((IDirectoryLocation)obj).FullName, Illumina.SecondaryAnalysis.Utilities.IsThisMono() ? StringComparison.Ordinal : StringComparison.OrdinalIgnoreCase);
            }
            else
            {
                return base.Equals(obj);
            }
        }

        public override int GetHashCode()
        {
            string fullName = Illumina.SecondaryAnalysis.Utilities.IsThisMono() ? FullName : FullName.ToUpper();
            return fullName.GetHashCode();
        }
    }
}
