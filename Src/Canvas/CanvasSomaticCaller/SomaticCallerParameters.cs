using System.Runtime.Serialization;
using System.Runtime.Serialization.Json;
using System.Security.Policy;

namespace somaticCallerParameters
{
    public class SomaticCallerParameters
    {
        [DataMember(Name = "MinimumVariantFrequenciesForInformativeSegment", EmitDefaultValue = true, IsRequired = true)]
        public int MinimumVariantFrequenciesForInformativeSegment { get; set; } = 50;
        [DataMember(Name = "MinimumCallSize ", EmitDefaultValue = true, IsRequired = true)]
        public int MinimumCallSize { get; set; } = 10;
        [DataMember(Name = "MaximumCopyNumber", EmitDefaultValue = true, IsRequired = true)]
        public int MaximumCopyNumber { get; set; } = 50000;       
        [DataMember(Name = "DefaultDeviationFactor", EmitDefaultValue = true, IsRequired = true)]
        public static float DefaultDeviationFactor { get; set; } = 1.75f;
        [DataMember(Name = "DefaultDeviationIndexCutoff", EmitDefaultValue = true, IsRequired = true)]
        public static int DefaultDeviationIndexCutoff { get; set; } = 18;
        [DataMember(Name = "DefaultCoverageWeighting", EmitDefaultValue = true, IsRequired = true)]
        public static double DefaultCoverageWeighting { get; set; } = 0.4;
        [DataMember(Name = "PrecisionWeightingFactor", EmitDefaultValue = true, IsRequired = true)]
        public double PrecisionWeightingFactor { get; set; } = 0.3333333333f;
        [DataMember(Name = "MaximumRelatedModels", EmitDefaultValue = true, IsRequired = true)]
        public static int MaximumRelatedModels { get; set; } = 5;
        [DataMember(Name = "PercentNormal2WeightingFactor", EmitDefaultValue = true, IsRequired = true)]
        public double PercentNormal2WeightingFactor { get; set; } = 0.275;
        [DataMember(Name = "DeviationScoreWeightingFactor", EmitDefaultValue = true, IsRequired = true)]
        public double DeviationScoreWeightingFactor { get; set; } = 0.275;
        [DataMember(Name = "CN2WeightingFactor", EmitDefaultValue = true, IsRequired = true)]
        public double CN2WeightingFactor { get; set; } = 0.375;     
        [DataMember(Name = "DiploidDistanceScoreWeightingFactor", EmitDefaultValue = true, IsRequired = true)]
        public double DiploidDistanceScoreWeightingFactor { get; set; } = 0.225;
        [DataMember(Name = "HeterogeneityScoreWeightingFactor", EmitDefaultValue = true, IsRequired = true)]
        public double HeterogeneityScoreWeightingFactor { get; set; } = 0.275;
        [DataMember(Name = "DeviationFactor", EmitDefaultValue = true, IsRequired = true)]
        public float DeviationFactor { get; set; } = 1.75f;
        [DataMember(Name = "DeviationIndexCutoff", EmitDefaultValue = true, IsRequired = true)]
        public int DeviationIndexCutoff { get; set; } = 18;
        [DataMember(Name = "CoverageWeighting", EmitDefaultValue = true, IsRequired = true)]
        public double CoverageWeighting { get; set; } = 0.4;
        [DataMember(Name = "MinAllowedPloidy", EmitDefaultValue = true, IsRequired = true)]
        public float MinAllowedPloidy { get; set; } = 0.5f;
        [DataMember(Name = "MaxAllowedPloidy", EmitDefaultValue = true, IsRequired = true)]
        public float MaxAllowedPloidy { get; set; } = 5.0f;
        [DataMember(Name = "UpperCoverageLevelWeightingFactor", EmitDefaultValue = true, IsRequired = true)]
        public float UpperCoverageLevelWeightingFactor { get; set; } = 2.5f;
        [DataMember(Name = "LowerCoverageLevelWeightingFactor", EmitDefaultValue = true, IsRequired = true)]
        public float LowerCoverageLevelWeightingFactor { get; set; } = 2.5f;
        [DataMember(Name = "CoverageLevelWeightingFactorStep", EmitDefaultValue = true, IsRequired = true)]
        public int CoverageLevelWeightingFactorLevels{ get; set; } = 80;
        [DataMember(Name = "UpperCentroidCutoff", EmitDefaultValue = true, IsRequired = true)]
        public float UpperCentroidCutoff { get; set; } = 2.5f;
        [DataMember(Name = "LowerCentroidCutoff", EmitDefaultValue = true, IsRequired = true)]
        public float LowerCentroidCutoff { get; set; } = 2.5f;
        [DataMember(Name = "CentroidCutoffStep", EmitDefaultValue = true, IsRequired = true)]
        public int CentroidCutoffStep { get; set; } = 80;
        [DataMember(Name = "DefaultCentroidCutoff", EmitDefaultValue = true, IsRequired = true)]
        public double DefaultCentroidCutoff { get; set; } = 0.03;
    }
}