using Isas.Shared.Utilities.FileSystem;

namespace CanvasNormalize
{
    interface IReferenceGenerator
    {
        void Run(IFileLocation outputFile);
    }
}
