using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Isas.Shared;

namespace CanvasNormalize
{
    interface IRatioCalculator
    {
        void Run(IFileLocation sampleBedFile, IFileLocation referenceBedFile, IFileLocation referencePloidyBedFile, IFileLocation outputBedFile);
    }
}
