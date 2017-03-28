using Illumina.Common.FileSystem;

namespace CanvasNormalize
{
    interface IReferenceGenerator
    {
        void Run(IFileLocation outputFile);
    }
}
