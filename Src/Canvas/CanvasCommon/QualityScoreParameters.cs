using System.Runtime.Serialization;
using System.Runtime.Serialization.Json;
using System.Security.Policy;

namespace CanvasCommon
{
    public class QualityScoreParameters
    {
        public double LogisticGermlineIntercept { get; set; } = -5.0123;
        public double LogisticGermlineLogBinCount { get; set; } = 4.9801;
        public double LogisticGermlineModelDistance { get; set; } = -5.5472;
        public double LogisticGermlineDistanceRatio { get; set; } = -1.7914;
        public double LogisticIntercept { get; set; } = -0.5143;
        public double LogisticLogBinCount { get; set; } = 0.8596;
        public double LogisticModelDistance { get; set; } = -50.4366;
        public double LogisticDistanceRatio { get; set; } = -0.6511;
        public double GeneralizedLinearFitIntercept { get; set; } = -3.65;
        public double GeneralizedLinearFitLogBinCount { get; set; } = -1.12;
        public double GeneralizedLinearFitModelDistance { get; set; } = 3.89;
        public double GeneralizedLinearFitMajorChromosomeCount { get; set; } = 0.47;
        public double GeneralizedLinearFitMafMean { get; set; } = -0.68;
        public double GeneralizedLinearFitLogMafCv { get; set; } = -0.25;
    }
}
