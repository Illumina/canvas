using Isas.Shared.Utilities.FileSystem;

namespace CanvasNormalize
{
    class CanvasNormalize
    {
        private readonly string CndFileSuffix = ".cnd"; // cnd: copy-number data
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
            CanvasNormalizeUtilities.WriteCndFile(parameters.tumorBedFile, parameters.weightedAverageNormalBedFile,
                ratios, new FileLocation(parameters.outBedFile.FullName + CndFileSuffix));
            
            return 0;
        }
    }
}
