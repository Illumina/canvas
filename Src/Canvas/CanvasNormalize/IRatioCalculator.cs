using System.Collections.Generic;
using CanvasCommon;
using Illumina.Common.FileSystem;

namespace CanvasNormalize
{
    interface IRatioCalculator
    {
        IEnumerable<SampleGenomicBin> Run(IFileLocation sampleBedFile, IFileLocation referenceBedFile);
    }
}
