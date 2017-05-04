namespace CanvasSomaticCaller
{
    public class SomaticCallerParameters
    {
        public int MinimumVariantFrequenciesForInformativeSegment { get; set; } = 50;
        public int MinimumCallSize { get; set; } = 10;
        public int MaximumCopyNumber { get; set; } = 50000;       
        public static float DefaultDeviationFactor { get; set; } = 1.75f;
        public static int DefaultDeviationIndexCutoff { get; set; } = 18;
        public static double DefaultCoverageWeighting { get; set; } = 0.4;
        public double PrecisionWeightingFactor { get; set; } = 0.3333333333f;
        public static int MaximumRelatedModels { get; set; } = 5;
        public double PercentNormal2WeightingFactor { get; set; } = 0.275;
        public double DeviationScoreWeightingFactor { get; set; } = 0.275;
        public double CN2WeightingFactor { get; set; } = 0.375;     
        public double DiploidDistanceScoreWeightingFactor { get; set; } = 0.225;
        public double HeterogeneityScoreWeightingFactor { get; set; } = 0.275;
        public float DeviationFactor { get; set; } = 1.75f;
        public int DeviationIndexCutoff { get; set; } = 18;
        public double CoverageWeighting { get; set; } = 0.4;
        public float MinAllowedPloidy { get; set; } = 0.5f;
        public float MaxAllowedPloidy { get; set; } = 5.0f;
        public float UpperCoverageLevelWeightingFactor { get; set; } = 2.5f;
        public float LowerCoverageLevelWeightingFactor { get; set; } = 2.5f;
        public int CoverageLevelWeightingFactorLevels{ get; set; } = 80;
        public float UpperCentroidCutoff { get; set; } = 2.5f;
        public float LowerCentroidCutoff { get; set; } = 2.5f;
        public int CentroidCutoffStep { get; set; } = 80;
        public double DefaultCentroidCutoff { get; set; } = 0.03;
        public double HeterogeneousClusterMedianCutoff { get; set; } = 1.25;
        public int HeterogeneousClustersCutoff { get; set; } = 1;
        public double DistanceRatio { get; set; } = 0.3;
        public double ClonalityIntercept { get; set; } = 2.774;
        public double ClonalityBestModelDistance { get; set; } = -13.575;
        public double ClonalityClusterEntropy { get; set; } = -3.522;
        public double ClonalityClusterMedianDistance { get; set; } = 2.808;
        public double ClonalityClusterMeanDistance { get; set; } = 6.8904;
        public double ClonalityClusterVariance { get; set; } = 14.9372;
        public double NumClusters { get; set; } = -0.0988;
        public double ModelDeviation { get; set; } = -7.3681;
        public double EvennessScoreThreshold { get; set; } = 94.5;
        public double CoverageWeightingWithMafSegmentation { get; set; } = 0.20;
    }
}