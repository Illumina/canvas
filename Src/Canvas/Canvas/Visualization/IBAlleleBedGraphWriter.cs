using Illumina.Common.FileSystem;
using Isas.Framework.DataTypes;

namespace Canvas.Visualization
{
    public interface IBAlleleBedGraphWriter
    {
        void Write(IFileLocation bafFile, BgzfFile bAllelesFile);
    }
}