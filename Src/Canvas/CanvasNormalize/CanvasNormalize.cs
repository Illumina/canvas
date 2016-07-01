using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using CanvasCommon;
using SequencingFiles;

namespace CanvasNormalize
{
    class CanvasNormalize
    {
        private IReferenceGenerator _referenceGenerator;
        private IRatioCalculator _ratioCalculator;
        public CanvasNormalize(IReferenceGenerator referenceGenerator, IRatioCalculator ratioCalculator)
        {
            _referenceGenerator = referenceGenerator;
            _ratioCalculator = ratioCalculator;
        }

        public int Run(CanvasNormalizeParameters parameters) 
        {
            _referenceGenerator.Run(parameters.weightedAverageNormalBedFile);
            var ratios = _ratioCalculator.Run(parameters.tumorBedFile, parameters.weightedAverageNormalBedFile);
            CanvasNormalizeUtilities.RatiosToCounts(ratios, parameters.ploidyBedFile, parameters.outBedFile);
            return 0;
        }
    }
}
