using Isas.Shared;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Illumina.SecondaryAnalysis
{
    public class TemporaryDirectory : DirectoryLocation, IDisposable
    {
        public TemporaryDirectory(string dir) : base(dir)
        {
            this.CreateClean();
        }

        public void Dispose()
        {
            if (!IsasConfiguration.GetConfiguration().RetainTempFiles)
                this.Delete();
        }
    }

    public static class TemporaryDirectoryLocationExtensions
    {
        public static TemporaryDirectory GetTemporaryDirectory(this IDirectoryLocation dirLocation, string newDirectory)
        {
            return new TemporaryDirectory(dirLocation.CreateSubdirectory(newDirectory).FullName);
        }
    }
}
