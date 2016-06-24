using System.Runtime.Serialization;
using System.Runtime.Serialization.Json;
using System.Security.Policy;

namespace CanvasCommon
{
    public class QualityScoreParameters
    {
        [DataMember(Name = "LogisticGermlineIntercept", EmitDefaultValue = true, IsRequired = true)]
        public double LogisticGermlineIntercept { get; set; } = -5.0123;
        [DataMember(Name = "LogisticGermlineLogBinCount", EmitDefaultValue = true, IsRequired = true)]
        public double LogisticGermlineLogBinCount { get; set; } = 4.9801;
        [DataMember(Name = "LogisticGermlineModelDistance", EmitDefaultValue = true, IsRequired = true)]
        public double LogisticGermlineModelDistance { get; set; } = -5.5472;
        [DataMember(Name = "LogisticGermlineDistanceRatio", EmitDefaultValue = true, IsRequired = true)]
        public double LogisticGermlineDistanceRatio { get; set; } = -1.7914;

        [DataMember(Name = "LogisticIntercept", EmitDefaultValue = true, IsRequired = true)]
        public double LogisticIntercept { get; set; } = -0.5143;
        [DataMember(Name = "LogisticLogBinCount", EmitDefaultValue = true, IsRequired = true)]
        public double LogisticLogBinCount { get; set; } = 0.8596;
        [DataMember(Name = "LogisticModelDistance", EmitDefaultValue = true, IsRequired = true)]
        public double LogisticModelDistance { get; set; } = -50.4366;
        [DataMember(Name = "LogisticDistanceRatio", EmitDefaultValue = true, IsRequired = true)]
        public double LogisticDistanceRatio { get; set; } = -0.6511;

        [DataMember(Name = "GeneralizedLinearFitIntercept", EmitDefaultValue = true, IsRequired = true)]
        public double GeneralizedLinearFitIntercept { get; set; } = -3.65;
        [DataMember(Name = "GeneralizedLinearFitLogBinCount", EmitDefaultValue = true, IsRequired = true)]
        public double GeneralizedLinearFitLogBinCount { get; set; } = -1.12;
        [DataMember(Name = "GeneralizedLinearFitModelDistance", EmitDefaultValue = true, IsRequired = true)]
        public double GeneralizedLinearFitModelDistance { get; set; } = 3.89;
        [DataMember(Name = "GeneralizedLinearFitMajorChromosomeCount", EmitDefaultValue = true, IsRequired = true)]
        public double GeneralizedLinearFitMajorChromosomeCount { get; set; } = 0.47;
        [DataMember(Name = "GeneralizedLinearFitMafMean", EmitDefaultValue = true, IsRequired = true)]
        public double GeneralizedLinearFitMafMean { get; set; } = -0.68;
        [DataMember(Name = "GeneralizedLinearFitLogMafCv", EmitDefaultValue = true, IsRequired = true)]
        public double GeneralizedLinearFitLogMafCv { get; set; } = -0.25;

    }
}
