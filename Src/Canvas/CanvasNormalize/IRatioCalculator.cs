using System.Collections.Generic;
using CanvasCommon;
using Isas.Shared.Utilities.FileSystem;

namespace CanvasNormalize
{
    interface IRatioCalculator
    {
        IEnumerable<GenomicBin> Run(IFileLocation sampleBedFile, IFileLocation referenceBedFile);
    }
}
