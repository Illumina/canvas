using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CanvasPedigreeCaller
{
    class CnvDistribution
    {
        private Array _probability;

        CnvDistribution(int nSamples, int nCopies, List<string> names)
        {
            if (names.Count != nSamples)
                throw new ArgumentException($"List of sample names must be equal to the number of samples ");
            var dimensionSizes = Enumerable.Repeat(nCopies, nSamples).ToArray();
            _probability = Array.CreateInstance(typeof(Int32), dimensionSizes);
        }

        public MarginalDistibution(string sample) { }
        public ConditionalDistibution(string sample) { }


    }
}
