using System;
using System.Linq;
using CanvasCommon;
using Isas.Manifests.NexteraManifest;

namespace CanvasNormalize
{
    internal class CanvasNormalizeFactory
    {
        private CanvasNormalizeParameters _parameters;
        private NexteraManifest _manifest;
        public CanvasNormalizeFactory(CanvasNormalizeParameters parameters)
        {
            _parameters = parameters;
            _manifest = _parameters.manifestFile == null ? null 
                : new NexteraManifest(_parameters.manifestFile.FullName, null, Console.WriteLine);
        }

        public IReferenceGenerator GetReferenceGenerator()
        {
            switch (_parameters.normalizationMode)
            {
                case CanvasNormalizeMode.BestLR2:
                    return new BestLR2ReferenceGenerator(_parameters.tumorBedFile, _parameters.normalBedFiles, _manifest);
                case CanvasNormalizeMode.WeightedAverage:
                    return new WeightedAverageReferenceGenerator(_parameters.normalBedFiles, _manifest);
                case CanvasNormalizeMode.PCA:
                    return new PCAReferenceGenerator(_parameters.tumorBedFile, _parameters.normalBedFiles.First(),
                        _manifest, _parameters.referenceBinCountRange.Min(), _parameters.referenceBinCountRange.Max());
                default:
                    throw new Exception(string.Format("Invalid CanvasNormalize mode '{0}'", _parameters.normalizationMode));
            }
        }

        public IRatioCalculator GetRatioCalculator()
        {
            switch (_parameters.normalizationMode)
            {
                case CanvasNormalizeMode.BestLR2:
                case CanvasNormalizeMode.WeightedAverage:
                    return new LSNormRatioCalculator(_manifest);
                case CanvasNormalizeMode.PCA:
                    return new RawRatioCalculator(_manifest, _parameters.referenceBinCountRange.Min(),
                        _parameters.referenceBinCountRange.Max());
                default:
                    throw new Exception(string.Format("Invalid CanvasNormalize mode '{0}'", _parameters.normalizationMode));
            }
        }
    }
}
