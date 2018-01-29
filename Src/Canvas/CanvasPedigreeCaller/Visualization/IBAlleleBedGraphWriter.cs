using Illumina.Common.FileSystem;
using Isas.Framework.DataTypes;

namespace CanvasPedigreeCaller.Visualization
{
    public interface IBAlleleBedGraphWriter
    {
        void Write(IFileLocation bafFile, BgzfFile bAllelesFile);
    }
}