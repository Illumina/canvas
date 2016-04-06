using System;
using System.Runtime.InteropServices;
using System.Text;
using Illumina.Common;
using ILMNcommon.Common;
using Microsoft.Win32.SafeHandles;
using Mono.Unix;
using Mono.Unix.Native;

namespace Isas.Shared.FileSystem
{
    public static class CrossPlatformLink
    {
        [DllImport("Kernel32.dll", CharSet = CharSet.Unicode, SetLastError = true)]
        private static extern bool CreateHardLink(string lpFileName, string lpExistingFileName, IntPtr lpSecurityAttributes);

        public static void HardlinkWindows(string source, string linkLocation)
        {
            if (!CreateHardLink(linkLocation, source, IntPtr.Zero))
                throw new ApplicationException(
                    $"Failed to create hardlink at '{linkLocation}' with source at '{source}'", CrossPlatform.GetLastWin32Exception());
        }


        [DllImport("kernel32.dll", CharSet = CharSet.Unicode, SetLastError = true)]
        private static extern bool CreateSymbolicLink(string lpSymlinkFileName, string lpTargetFileName, int dwFlags);

        public static void SymlinkWindows(string source, string linkLocation)
        {
            if (!CreateSymbolicLink(linkLocation, source, 0))
                throw new ApplicationException(
                    $"Failed to create symbolic link at '{linkLocation}' pointing to '{source}'", CrossPlatform.GetLastWin32Exception());
        }

        public static void HardlinkUnix(string source, string linkLocation)
        {
            int r = Syscall.link(source, linkLocation);
            // wrap any exception with a more friendly error message
            Action a = () => UnixMarshal.ThrowExceptionForLastErrorIf(r);
            a.WrapAnyExceptionWith(
               nestedException => new ApplicationException(
                   $"Failed to create hardlink at '{linkLocation}' with source at '{source}'", nestedException));
        }

        public static void SymlinkUnix(string source, string linkLocation)
        {
            int r = Syscall.symlink(source, linkLocation);
            // wrap any exception with a more friendly error message
            Action a = () => UnixMarshal.ThrowExceptionForLastErrorIf(r);
            a.WrapAnyExceptionWith(
                nestedException => new ApplicationException(
                    $"Failed to create symbolic link at '{linkLocation}' pointing to '{source}'", nestedException));
        }

        public static void Symlink(string source, string linkLocation)
        {
            if (Utilities.IsThisMono())
                SymlinkUnix(source, linkLocation);
            else
                SymlinkWindows(source, linkLocation);
        }

        public static void Hardlink(string source, string linkLocation)
        {
            if (Utilities.IsThisMono())
                HardlinkUnix(source, linkLocation);
            else
                HardlinkWindows(source, linkLocation);
        }

        public static string ReadLink(string path)
        {
            if (Utilities.IsThisMono())
                return ReadLinkUnix(path);
            else
                return ReadLinkWindows(path);
        }

        /// <summary>
        /// This is placed in a separate method call so that the CLR does not attempt to load the Mono.Posix assembly on Windows.
        /// We don't include the Mono.Posix assembly as part of the build. It is provided automatically When running under mono
        /// </summary>
        private static string ReadLinkUnix(string path)
        {
            try
            {
                return UnixPath.GetCompleteRealPath(path);
            }
            catch (Exception e)
            {
                throw new ApplicationException($"Unable to readlink path '{path}'", e);
            }
        }

        private const int FILE_SHARE_WRITE = 2;
        private const int CREATION_DISPOSITION_OPEN_EXISTING = 3;
        private const int FILE_FLAG_BACKUP_SEMANTICS = 0x02000000;
        [DllImport("kernel32.dll", CharSet = CharSet.Unicode, SetLastError = true)]
        private static extern int GetFinalPathNameByHandle(IntPtr handle, [In, Out] StringBuilder path, int bufLen, int flags);
        [DllImport("kernel32.dll", EntryPoint = "CreateFileW", CharSet = CharSet.Unicode, SetLastError = true)]
        private static extern SafeFileHandle CreateFile(string lpFileName, int dwDesiredAccess, int dwShareMode,
 IntPtr securityAttributes, int dwCreationDisposition, int dwFlagsAndAttributes, IntPtr hTemplateFile);
        private static string ReadLinkWindows(string path)
        {
            int maxReadLinkPathCapacity = 32767;
            var readLinkPath = new StringBuilder(512, maxReadLinkPathCapacity);
            int returnedLength;
            using (
                SafeFileHandle handle = CreateFile(path, 0, FILE_SHARE_WRITE, IntPtr.Zero, CREATION_DISPOSITION_OPEN_EXISTING, FILE_FLAG_BACKUP_SEMANTICS,
                    IntPtr.Zero))
            {
                if (handle == null || handle.IsInvalid)
                    throw new InvalidOperationException($"Could not get handle to path at '{path}'", CrossPlatform.GetLastWin32Exception());
                returnedLength = GetFinalPathNameByHandle(handle.DangerousGetHandle(), readLinkPath, maxReadLinkPathCapacity, 0);
            }
            if (returnedLength <= 0 || returnedLength - 1 > maxReadLinkPathCapacity)
                throw new InvalidOperationException($"Readlink error. Could not get full path for '{path}'", CrossPlatform.GetLastWin32Exception());

            // this method returns some weird prefix that we need to fix:
            // network path: \\?\UNC\sd-isilon\bioinfoSD --> \\sd-isilon\bioinfoSD 
            // local path: \\?\C:\Projects --> C:\Projects
            string readLinkResult = readLinkPath.ToString();
            string localPrefix = @"\\?\";
            string networkPrefix = localPrefix + @"UNC\";
            if (readLinkResult.StartsWith(networkPrefix))
                return $@"\\{readLinkResult.Substring(networkPrefix.Length)}";
            if (readLinkResult.StartsWith(localPrefix))
                return readLinkResult.Substring(localPrefix.Length);
            return readLinkResult;
        }

        public static string GetRelativePath(string startingPath, string endingPath)
        {
            var startingUri = new Uri(startingPath);
            var endingUri = new Uri(endingPath);
            var relativeUri = startingUri.MakeRelativeUri(endingUri);
            return relativeUri.ToString();
        }
    }
}
