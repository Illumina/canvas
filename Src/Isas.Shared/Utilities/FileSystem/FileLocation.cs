using System;
using System.IO;
using System.Security.AccessControl;
using ILMNcommon.Common;
using Isas.Shared.FileSystem;

namespace Isas.Shared
{
    /// <summary>
    /// Provides an interface for properties and instance methods for the creation, copying, deletion, moving, and opening of files, and aids in the creation of FileStream objects. This class cannot be inherited.
    /// </summary>
    public interface IFileLocation
    {
        /// <summary>
        /// Permanently deletes a file.
        /// </summary>
        void Delete();

        /// <summary>
        /// Refreshes the state of the object.
        /// </summary>
        void Refresh();

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
        /// The creation date and time of the current FileSystemInfo object.
        /// </value>
        DateTime CreationTime { get; set; }

        /// <summary>
        /// Gets or sets the creation time, in coordinated universal time (UTC), of the current file or directory.
        /// </summary>
        /// <value>
        /// The creation date and time in UTC format of the current FileSystemInfo object.
        /// </value>
        DateTime CreationTimeUtc { get; set; }

        /// <summary>
        /// Gets a value indicating whether a file exists.
        /// </summary>
        /// <value>
        ///   <c>true</c> if the file exists; <c>false</c> if the file does not exist or if the file is a directory.
        /// </value>
        bool Exists { get; }

        /// <summary>
        /// Gets the string representing the extension part of the file.
        /// </summary>
        /// <value>
        /// A string containing the FileSystemInfo extension.
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
        /// Gets the full path of the file after following all symlinks and removing any intermediate relative paths (e.g. "/../" or "/./")
        /// </summary>
        /// <value>
        /// A string containing the full path.
        /// </value>
        string FullNameCanonical { get; }

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
        /// Gets the name of the file without the directory
        /// </summary>
        /// <value>
        /// The name of the file.
        /// </value>
        string Name { get; }

        /// <summary>
        /// Gets an instance of the parent directory
        /// </summary>
        /// <value>
        /// An IDirectoryLocation object representing the parent directory of this file.
        /// </value>
        IDirectoryLocation Directory { get; }

        /// <summary>
        /// Gets a string representing the directory's full path.
        /// </summary>
        /// <value>
        /// A string representing the directory's full path.
        /// </value>
        string DirectoryName { get; }

        /// <summary>
        /// Gets or sets a value that determines if the current file is read only.
        /// </summary>
        /// <value>
        /// <c>true</c> if the current file is read only; otherwise, <c>false</c>.
        /// </value>
        bool IsReadOnly { get; set; }

        /// <summary>
        /// Gets the size, in bytes, of the current file.
        /// </summary>
        /// <value>
        /// The size of the current file in bytes.
        /// </value>
        long Length { get; }

        /// <summary>
        /// Creates a StreamWriter that appends text to the file represented by this instance of the IFileLocation.
        /// </summary>
        /// <returns>A new StreamWriter.</returns>
        StreamWriter AppendText();

        /// <summary>
        /// Copies an existing file to a new file, overwriting any existing file.
        /// </summary>
        /// <param name="destination">The location of the new file to copy to.</param>
        /// <returns>
        /// The location of the new file
        /// </returns>
        IFileLocation CopyTo(IFileLocation destination);

        /// <summary>
        /// Copies an existing file into a parent directory, overwriting any existing file.
        /// </summary>
        /// <param name="destination">The parent directory of the new file to copy to.</param>
        /// <returns>
        /// The location of the new file
        /// </returns>
        IFileLocation CopyInto(IDirectoryLocation destination);

        /// <summary>
        /// Moves an existing file into a parent directory, overwriting any existing file.
        /// </summary>
        /// <param name="destination">The parent directory for the new location of the file.</param>
        /// <returns>
        /// The new location of the file
        /// </returns>
        IFileLocation MoveInto(IDirectoryLocation destination);

        /// <summary>
        /// Moves an existing file to a new location, overwriting any existing file.
        /// If the new location is on a different volume, this is equivalent to CopyTo + Delete
        /// </summary>
        /// <param name="destination">The new location of the file.</param>
        /// <returns>
        /// The new location of the file
        /// </returns>
        IFileLocation MoveTo(IFileLocation destination);

        /// <summary>
        /// Creates a file.
        /// </summary>
        /// <returns>A new file.</returns>
        Stream Create();

        /// <summary>
        /// Creates a StreamWriter that writes a new text file.
        /// </summary>
        /// <returns>A new StreamWriter.</returns>
        StreamWriter CreateText();

        /// <summary>
        /// Decrypts a file that was encrypted by the current account using the Encrypt method.
        /// </summary>
        void Decrypt();

        /// <summary>
        /// Encrypts a file so that only the account used to encrypt the file can decrypt it.
        /// </summary>
        void Encrypt();

        /// <summary>
        /// Gets a FileSecurity object that encapsulates the access control list (ACL) entries for the file described by the current IFileLocation object.
        /// </summary>
        /// <returns>A FileSecurity object that encapsulates the access control rules for the current file.</returns>
        FileSecurity GetAccessControl();

        /// <summary>
        /// Gets a FileSecurity object that encapsulates the specified type of access control list (ACL) entries for the file described by the current IFileLocation object.
        /// </summary>
        /// <param name="includeSections">One of the AccessControlSections values that specifies which group of access control entries to retrieve.</param>
        /// <returns>A FileSecurity object that encapsulates the access control rules for the current file. </returns>
        FileSecurity GetAccessControl(AccessControlSections includeSections);

        /// <summary>
        /// Opens a file in the specified mode.
        /// </summary>
        /// <param name="mode">A FileMode constant specifying the mode (for example, Open or Append) in which to open the file.</param>
        /// <returns>A file opened in the specified mode, with read/write access and unshared.</returns>
        Stream Open(FileMode mode);

        /// <summary>
        /// Opens a file in the specified mode with read, write, or read/write access.
        /// </summary>
        /// <param name="mode">A FileMode constant specifying the mode (for example, Open or Append) in which to open the file.</param>
        /// <param name="access">A FileAccess constant specifying whether to open the file with Read, Write, or ReadWrite file access.</param>
        /// <returns>A FileStream object opened in the specified mode and access, and unshared.</returns>
        Stream Open(FileMode mode, FileAccess access);

        /// <summary>
        /// Opens a file in the specified mode with read, write, or read/write access and the specified sharing option.
        /// </summary>
        /// <param name="mode">A FileMode constant specifying the mode (for example, Open or Append) in which to open the file.</param>
        /// <param name="access">A FileAccess constant specifying whether to open the file with Read, Write, or ReadWrite file access.</param>
        /// <param name="share">A FileShare constant specifying the type of access other FileStream objects have to this file.</param>
        /// <returns>A FileStream object opened with the specified mode, access, and sharing options.</returns>
        Stream Open(FileMode mode, FileAccess access, FileShare share);

        /// <summary>
        /// Creates a read-only FileStream.
        /// </summary>
        /// <returns>A new read-only FileStream object.</returns>
        Stream OpenRead();

        /// <summary>
        /// Creates a StreamReader with UTF8 encoding that reads from an existing text file.
        /// </summary>
        /// <returns>A new StreamReader with UTF8 encoding.</returns>
        StreamReader OpenText();

        /// <summary>
        /// Creates a write-only FileStream.
        /// </summary>
        /// <returns>A write-only unshared FileStream object for a new or existing file.</returns>
        Stream OpenWrite();

        /// <summary>
        /// Applies access control list (ACL) entries described by a FileSecurity object to the file described by the current IFileLocation object.
        /// </summary>
        /// <param name="fileSecurity">A FileSecurity object that describes an access control list (ACL) entry to apply to the current file.</param>
        void SetAccessControl(FileSecurity fileSecurity);

        IFileLocation CreateAbsoluteSymlinkAt(IFileLocation linkLocation);
        bool TryCreateAbsoluteSymlinkAt(IFileLocation linkLocation, out Exception e);
        IFileLocation CreateAbsoluteSymlinkOrCopy(IFileLocation linkLocation);
        IFileLocation CreateHardlinkAt(IFileLocation linkLocation);
        bool TryCreateHardlinkAt(IFileLocation linkLocation, out Exception e);
        IFileLocation CreateHardlinkOrCopy(IFileLocation linkLocation);
        IFileLocation CreateRelativeSymlinkAt(IFileLocation linkLocation);
        bool TryCreateRelativeSymlinkAt(IFileLocation linkLocation, out Exception e);
        IFileLocation CreateRelativeSymlinkOrCopy(IFileLocation linkLocation);
    }

    public static class FileLocationExtensions
    // Additional extensions are found in the SecondaryAnalysis.Utilities class
    {
        /// <summary>
        /// Ensures an instance of the file in the file system exists. Like Linux "touch".
        /// </summary>
        /// <param name="file"></param>
        public static IFileLocation Touch(this IFileLocation file)
        {
            if (!file.Exists)
            {
                if (!file.Directory.Exists)
                    file.Directory.Create();

                file.Create().Dispose();
            }
            file.Refresh();
            return file;
        }

        /// <summary>
        /// Returns a new instance of a file with the <paramref name="name"/> appended to the end.
        /// </summary>
        /// <param name="file"></param>
        /// <param name="name"></param>
        /// <returns></returns>
        public static IFileLocation AppendName(this IFileLocation file, string name)
        {
            return file.Directory.GetFileLocation(file.Name + name);
        }

        /// <summary>
        /// Not making this an extension method so that it doesn't collide with the instance method
        /// We still want to be able to refer to it statically though
        /// </summary>
        /// <param name="source"></param>
        /// <param name="destination"></param>
        public static void MoveTo(IFileLocation source, IFileLocation destination)
        {
            source.MoveTo(destination);
        }

        /// <summary>
        /// void return type
        /// </summary>
        public static void MoveAndLinkAction(this IFileLocation source, IFileLocation destination)
        {
            source.MoveAndLink(destination);
        }

        /// <summary>
        /// Move a file to the destination, then put a symlink in the original location
        /// pointing to the destination copy.
        /// </summary>
        /// <param name="sourcePath"></param>
        /// <param name="destinationPath"></param>
        public static IFileLocation MoveAndLink(this IFileLocation source, IFileLocation destination)
        {
            // If the destination is null, return null
            if (destination == null)
                return null;

            source.MoveTo(destination);
            destination.CreateRelativeSymlinkOrCopy(source);

            return destination;
        }

        public static IFileLocation MoveAndLinkInto(this IFileLocation source, IDirectoryLocation destination)
        {
            return source.MoveAndLink(destination.GetFileLocation(source.Name));
        }

        /// <summary>
        /// A static delegate for IFileLocation.CreateHardlinkOrCopy
        /// </summary>
        public static void CreateHardlinkOrCopy(IFileLocation source, IFileLocation linkLocation)
        {
            source.CreateHardlinkOrCopy(linkLocation);
        }

        /// <summary>
        /// A static delegate for IFileLocation.CreateHardlinkAt
        /// </summary>
        public static void CreateHardlink(IFileLocation source, IFileLocation linkLocation)
        {
            source.CreateHardlinkAt(linkLocation);
        }

        /// <summary>
        /// Creates a relative symlink in <paramref name="linkParentDirectory"/> pointing to <paramref name="sourceFile"/>.
        /// Symlinks in results should be relative to allow browsing with different mount-point (e.g. windows)
        /// </summary>
        /// <param name="sourceFile">Source file</param>
        /// <param name="linkParentDirectory">Folder in which to create the relative symlink.</param>
        /// <returns>IFileLocation representing the new symlink.</returns>
        public static IFileLocation CreateRelativeSymlinkOrCopyIn(this IFileLocation sourceFile, IDirectoryLocation linkParentDirectory)
        {
            IFileLocation dest = linkParentDirectory.GetFileLocation(sourceFile.Name);
            return sourceFile.CreateRelativeSymlinkOrCopy(dest);
        }

        public static IFileLocation ReplaceEnd(this IFileLocation file, string oldSuffix, string newSuffix)
        {
            return file.Directory.GetFileLocation(file.Name.ReplaceEnd(oldSuffix, newSuffix));
        }

        public static IFileLocation CreateAbsoluteSymlinkOrCopyIn(this IFileLocation sourceFile, IDirectoryLocation destFolder)
        {
            IFileLocation dest = destFolder.GetFileLocation(sourceFile.Name);
            return sourceFile.CreateAbsoluteSymlinkOrCopy(dest);
        }
    }

    /// <summary>
    /// Provides properties and instance methods for the creation, copying, deletion, moving, and opening of files, and aids in the creation of FileStream objects.
    /// </summary>
    public class FileLocation : FileSystemLocationBase, IFileLocation
    {
        private readonly FileInfo _instance;

        /// <summary>
        /// Initializes a new instance of the <see cref="FileLocation"/> class.
        /// </summary>
        /// <param name="fileInfo">The FileInfo object to wrap.</param>
        public FileLocation(FileInfo fileInfo)
        {
            if (fileInfo == null)
                throw new ArgumentNullException(nameof(fileInfo));
            _instance = fileInfo;
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="FileLocation"/> class.
        /// </summary>
        /// <param name="filePath">The file path.</param>
        public FileLocation(string filePath) : this(new FileInfo(filePath))
        {
        }

        /// <summary>
        /// Permanently deletes a file.
        /// </summary>
        public override void Delete()
        {
            _instance.Delete();
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
            get { return _instance.Attributes; }
            set { _instance.Attributes = value; }
        }

        /// <summary>
        /// Gets or sets the creation time of the current file or directory.
        /// </summary>
        /// <value>
        /// The creation date and time of the current FileSystemInfo object.
        /// </value>
        public override DateTime CreationTime
        {
            get { return _instance.CreationTime; }
            set { _instance.CreationTime = value; }
        }

        /// <summary>
        /// Gets or sets the creation time, in coordinated universal time (UTC), of the current file or directory.
        /// </summary>
        /// <value>
        /// The creation date and time in UTC format of the current FileSystemInfo object.
        /// </value>
        public override DateTime CreationTimeUtc
        {
            get { return _instance.CreationTimeUtc; }
            set { _instance.CreationTimeUtc = value; }
        }

        /// <summary>
        /// Gets a value indicating whether a file exists.
        /// </summary>
        /// <value>
        /// <c>true</c> if the file exists; <c>false</c> if the file does not exist or if the file is a directory.
        /// </value>
        public override bool Exists
        {
            get { Refresh(); return _instance.Exists; }
        }

        /// <summary>
        /// Gets the string representing the extension part of the file.
        /// </summary>
        /// <value>
        /// A string containing the FileSystemInfo extension.
        /// </value>
        public override string Extension
        {
            get { return _instance.Extension; }
        }

        /// <summary>
        /// Gets the full path of the directory or file that was used to contruct this FileLocation.
        /// </summary>
        /// <value>
        /// A string containing the full path.
        /// </value>
        public override string FullName
        {
            get { return _instance.FullName; }
        }

        /// <summary>
        /// The real path to this file after following any symlinks and removing any intermediate relative paths (e.g. "/../" or "/./")
        /// The actual file does not need to exist
        /// </summary>
        public string FullNameCanonical
        {
            get
            {
                if (!Exists)
                    return Path.Combine(Directory.FullNameCanonical, Name);
                return CrossPlatformLink.ReadLink(FullName);
            }
        }

        /// <summary>
        /// Gets or sets the time the current file or directory was last accessed.
        /// </summary>
        /// <value>
        /// The time that the current file or directory was last accessed.
        /// </value>
        public override DateTime LastAccessTime
        {
            get { return _instance.LastAccessTime; }
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
            get { return _instance.LastAccessTimeUtc; }
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
            get { return _instance.LastWriteTime; }
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
            get { return _instance.LastWriteTimeUtc; }
            set { _instance.LastWriteTimeUtc = value; }
        }

        /// <summary>
        /// Gets the name of the file.
        /// </summary>
        /// <value>
        /// The name of the file.
        /// </value>
        public override string Name
        {
            get { return _instance.Name; }
        }

        /// <summary>
        /// Creates a StreamWriter that appends text to the file represented by this instance of the FileLocation.
        /// </summary>
        /// <returns>
        /// A new StreamWriter.
        /// </returns>
        public StreamWriter AppendText()
        {
            return _instance.AppendText();
        }

        /// <summary>
        /// Copies an existing file to a new file, overwriting any existing file.
        /// </summary>
        /// <param name="destination">The location of the new file to copy to.</param>
        /// <returns>
        /// The location of the new file
        /// </returns>
        public IFileLocation CopyTo(IFileLocation destination)
        {
            if (FullNameCanonical == destination.FullNameCanonical)
            {
                if (!Exists)
                    throw new ApplicationException($"Source and destination are the same, but source doesn't exist at {FullName}");
                return destination;
            }

            destination.Directory.Create();
            destination.Delete();
            File.Copy(FullName, destination.FullName);
            return destination;
        }

        /// <summary>
        /// Copies an existing file into a parent directory, overwriting any existing file.
        /// </summary>
        /// <param name="destination">The parent directory of the new file to copy to.</param>
        /// <returns>
        /// The location of the new file
        /// </returns>
        public IFileLocation CopyInto(IDirectoryLocation destination)
        {
            var result = destination.GetFileLocation(Name);
            CopyTo(result);
            return result;
        }

        /// <summary>
        /// Moves an existing file into a parent directory, overwriting any existing file.
        /// </summary>
        /// <param name="destination">The parent directory for the new location of the file.</param>
        /// <returns>
        /// The new location of the file
        /// </returns>
        public IFileLocation MoveInto(IDirectoryLocation destination)
        {
            var result = destination.GetFileLocation(Name);
            MoveTo(result);
            return result;
        }

        /// <summary>
        /// Moves an existing file to a new location, overwriting any existing file.
        /// If the new location is on a different volume, this is equivalent to CopyTo + Delete
        /// </summary>
        /// <param name="destination">The new location of the file.</param>
        /// <returns>
        /// The new location of the file
        /// </returns>
        public IFileLocation MoveTo(IFileLocation destination)
        {
            if (FullNameCanonical == destination.FullNameCanonical)
            {
                if (!Exists)
                    throw new ApplicationException($"Source and destination are the same, but source doesn't exist at {FullName}");
                return destination;
            }

            destination.Directory.Create();
            destination.Delete();
            File.Move(FullNameCanonical, destination.FullName);
            return destination;
        }

        /// <summary>
        /// Creates a file.
        /// </summary>
        /// <returns>
        /// A new file.
        /// </returns>
        public Stream Create()
        {
            return _instance.Create();
        }

        /// <summary>
        /// Creates a StreamWriter that writes a new text file.
        /// </summary>
        /// <returns>
        /// A new StreamWriter.
        /// </returns>
        public StreamWriter CreateText()
        {
            return _instance.CreateText();
        }

        /// <summary>
        /// Decrypts a file that was encrypted by the current account using the Encrypt method.
        /// </summary>
        public void Decrypt()
        {
            _instance.Decrypt();
        }

        /// <summary>
        /// Encrypts a file so that only the account used to encrypt the file can decrypt it.
        /// </summary>
        public void Encrypt()
        {
            _instance.Encrypt();
        }

        /// <summary>
        /// Gets a FileSecurity object that encapsulates the access control list (ACL) entries for the file described by the current FileLocation object.
        /// </summary>
        /// <returns>
        /// A FileSecurity object that encapsulates the access control rules for the current file.
        /// </returns>
        public FileSecurity GetAccessControl()
        {
            return _instance.GetAccessControl();
        }

        /// <summary>
        /// Gets a FileSecurity object that encapsulates the specified type of access control list (ACL) entries for the file described by the current FileLocation object.
        /// </summary>
        /// <param name="includeSections">One of the AccessControlSections values that specifies which group of access control entries to retrieve.</param>
        /// <returns>
        /// A FileSecurity object that encapsulates the access control rules for the current file.
        /// </returns>
        public FileSecurity GetAccessControl(AccessControlSections includeSections)
        {
            return _instance.GetAccessControl(includeSections);
        }

        /// <summary>
        /// Opens a file in the specified mode.
        /// </summary>
        /// <param name="mode">A FileMode constant specifying the mode (for example, Open or Append) in which to open the file.</param>
        /// <returns>
        /// A file opened in the specified mode, with read/write access and unshared.
        /// </returns>
        public Stream Open(FileMode mode)
        {
            return _instance.Open(mode);
        }

        /// <summary>
        /// Opens a file in the specified mode with read, write, or read/write access.
        /// </summary>
        /// <param name="mode">A FileMode constant specifying the mode (for example, Open or Append) in which to open the file.</param>
        /// <param name="access">A FileAccess constant specifying whether to open the file with Read, Write, or ReadWrite file access.</param>
        /// <returns>
        /// A FileStream object opened in the specified mode and access, and unshared.
        /// </returns>
        public Stream Open(FileMode mode, FileAccess access)
        {
            return _instance.Open(mode, access);
        }

        /// <summary>
        /// Opens a file in the specified mode with read, write, or read/write access and the specified sharing option.
        /// </summary>
        /// <param name="mode">A FileMode constant specifying the mode (for example, Open or Append) in which to open the file.</param>
        /// <param name="access">A FileAccess constant specifying whether to open the file with Read, Write, or ReadWrite file access.</param>
        /// <param name="share">A FileShare constant specifying the type of access other FileStream objects have to this file.</param>
        /// <returns>
        /// A FileStream object opened with the specified mode, access, and sharing options.
        /// </returns>
        public Stream Open(FileMode mode, FileAccess access, FileShare share)
        {
            return _instance.Open(mode, access, share);
        }

        /// <summary>
        /// Creates a read-only FileStream.
        /// </summary>
        /// <returns>
        /// A new read-only FileStream object.
        /// </returns>
        public Stream OpenRead()
        {
            return _instance.OpenRead();
        }

        /// <summary>
        /// Creates a StreamReader with UTF8 encoding that reads from an existing text file.
        /// </summary>
        /// <returns>
        /// A new StreamReader with UTF8 encoding.
        /// </returns>
        public StreamReader OpenText()
        {
            return _instance.OpenText();
        }

        /// <summary>
        /// Creates a write-only FileStream.
        /// </summary>
        /// <returns>
        /// A write-only unshared FileStream object for a new or existing file.
        /// </returns>
        public Stream OpenWrite()
        {
            return _instance.OpenWrite();
        }

        /// <summary>
        /// Applies access control list (ACL) entries described by a FileSecurity object to the file described by the current FileLocation object.
        /// </summary>
        /// <param name="fileSecurity">A FileSecurity object that describes an access control list (ACL) entry to apply to the current file.</param>
        public void SetAccessControl(FileSecurity fileSecurity)
        {
            _instance.SetAccessControl(fileSecurity);
        }

        /// <summary>
        /// Gets an instance of the parent directory
        /// </summary>
        /// <value>
        /// An IDirectoryLocation object representing the parent directory of this file.
        /// </value>
        public IDirectoryLocation Directory
        {
            get { return new DirectoryLocation(_instance.Directory); }
        }

        /// <summary>
        /// Gets a string representing the directory's full path.
        /// </summary>
        /// <value>
        /// A string representing the directory's full path.
        /// </value>
        public string DirectoryName
        {
            get { return _instance.DirectoryName; }
        }

        /// <summary>
        /// Gets or sets a value that determines if the current file is read only.
        /// </summary>
        /// <value>
        /// <c>true</c> if the current file is read only; otherwise, <c>false</c>.
        /// </value>
        public bool IsReadOnly
        {
            get { return _instance.IsReadOnly; }
            set { _instance.IsReadOnly = value; }
        }

        /// <summary>
        /// Gets the size, in bytes, of the current file.
        /// </summary>
        /// <value>
        /// The size of the current file in bytes.
        /// </value>
        public long Length
        {
            get { return _instance.Length; }
        }

        public override string ToString()
        {
            return FullName;
        }

        public override bool Equals(object obj)
        {
            var item = obj as IFileLocation;
            if (item == null) return false;
            return FullName.Equals(item.FullName, Utilities.IsThisMono() ? StringComparison.Ordinal : StringComparison.OrdinalIgnoreCase);
        }

        public override int GetHashCode()
        {
            string fullName = Utilities.IsThisMono() ? FullName : FullName.ToUpper();
            return fullName.GetHashCode();
        }

        public bool TryCreateAbsoluteSymlinkAt(IFileLocation linkLocation, out Exception e)
        {
            e = null;
            try
            {
                CreateAbsoluteSymlinkAt(linkLocation);
                return true;
            }
            catch (Exception ex)
            {
                e = ex;
                return false;
            }
        }

        public bool TryCreateRelativeSymlinkAt(IFileLocation linkLocation, out Exception e)
        {
            e = null;
            try
            {
                CreateRelativeSymlinkAt(linkLocation);
                return true;
            }
            catch (Exception ex)
            {
                e = ex;
                return false;
            }
        }

        public bool TryCreateHardlinkAt(IFileLocation linkLocation, out Exception e)
        {
            e = null;
            try
            {
                CreateHardlinkAt(linkLocation);
                return true;
            }
            catch (Exception ex)
            {
                e = ex;
                return false;
            }
        }

        public IFileLocation CreateAbsoluteSymlinkAt(IFileLocation linkLocation)
        {
            linkLocation.Delete();
            CrossPlatformLink.Symlink(FullName, linkLocation.FullName);
            return linkLocation;
        }

        public IFileLocation CreateRelativeSymlinkAt(IFileLocation linkLocation)
        {
            linkLocation.Delete();
            CrossPlatformLink.Symlink(CrossPlatformLink.GetRelativePath(linkLocation.FullName, FullName), linkLocation.FullName);
            return linkLocation;
        }

        public IFileLocation CreateHardlinkAt(IFileLocation linkLocation)
        {
            linkLocation.Directory.Create();
            linkLocation.Delete();
            CrossPlatformLink.Hardlink(FullName, linkLocation.FullName);
            return linkLocation;
        }

        public IFileLocation CreateHardlinkIn(IDirectoryLocation directory)
        {
            var linkLocation = directory.GetFileLocation(Name);
            return CreateHardlinkAt(linkLocation);
        }

        /// <summary>
        /// Symlinks an existing file to a new location.
        /// If the symlink cannot be created (Windows may require elevated privileges to do this),
        ///  then a copy is created instead if <paramref name="copyOnFail"/> is true.
        /// </summary>
        /// <param name="sourceFile">The instance of IFileLocation this method extends</param>
        /// <param name="linkLocation">The IFileLocation representing the symlink, 
        /// which points to <paramref name="sourceFile"/></param>
        /// <param name="copyOnFail">If true, copies the file if a symlink cannot be created.</param>
        /// <returns>An IFileLocation of the destination file, or null if the symlink failed and <paramref name="copyOnFail"/> is false.</returns>
        public IFileLocation CreateAbsoluteSymlinkOrCopy(IFileLocation linkLocation)
        {
            Action a = () => CreateAbsoluteSymlinkAt(linkLocation);
            if (!a.Try())
                CopyTo(linkLocation);
            return linkLocation;
        }

        public IFileLocation CreateHardlinkOrCopy(IFileLocation linkLocation)
        {
            Action a = () => CreateHardlinkAt(linkLocation);
            if (!a.Try())
                CopyTo(linkLocation);
            return linkLocation;
        }

        /// <summary>
        /// Creates a relative symlink at <paramref name="linkLocation"/> pointing to <paramref name="sourceFile"/>
        /// Symlinks in results should be relative to allow browsing with different mount-point (e.g. windows)
        /// </summary>
        /// <param name="sourceFile">Source file</param>
        /// <param name="linkLocation">Symbolic link to be created.</param>
        /// <returns>IFileLocation representing the new symlink.</returns>
        public IFileLocation CreateRelativeSymlinkOrCopy(IFileLocation linkLocation)
        {
            Action a = () => CreateRelativeSymlinkAt(linkLocation);
            if (!a.Try())
                CopyTo(linkLocation);
            return linkLocation;
        }
    }
}
