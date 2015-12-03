using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace Isas.Shared
{
    internal static class Converters
    {
        // Note - symlinks will look just like regular FileInfo/DirectoryInfo objects
        internal static FileSystemLocationBase[] WrapFileSystemInfos(this IEnumerable<FileSystemInfo> input)
        {
            return input
                .Select<FileSystemInfo, FileSystemLocationBase>(item =>
                {
                    if (item is FileInfo)
                        return new FileLocation(item as FileInfo);

                    if (item is DirectoryInfo)
                        return new DirectoryLocation(item as DirectoryInfo);

                    throw new NotImplementedException(string.Format(
                        CultureInfo.InvariantCulture,
                        "The type {0} is not recognized by the System.IO.Abstractions library.",
                        item.GetType().AssemblyQualifiedName
                    ));
                })
                .ToArray();
        }

        internal static IDirectoryLocation[] WrapDirectories(this IEnumerable<DirectoryInfo> input)
        {
            return input.Select(f => new DirectoryLocation(f)).ToArray();
        }

        internal static IFileLocation[] WrapFiles(this IEnumerable<FileInfo> input)
        {
            return input.Select(f => new FileLocation(f)).ToArray();
        }
    }
}
