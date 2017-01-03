using System.Collections.Generic;
using CanvasCommon;
using Isas.Shared.Utilities.FileSystem;

namespace CanvasNormalize
{
    interface IRatioCalculator
    {
        IEnumerable<SampleGenomicBin> Run(IFileLocation sampleBedFile, IFileLocation referenceBedFile);
    }
}
