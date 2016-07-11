using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Isas.Shared;
using CanvasCommon;

namespace CanvasNormalize
{
    interface IRatioCalculator
    {
        IEnumerable<GenomicBin> Run(IFileLocation sampleBedFile, IFileLocation referenceBedFile);
    }
}
