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
        IEnumerable<Tuple<string, int, int, double>> Run(IFileLocation sampleBedFile, IFileLocation referenceBedFile);
    }
}
