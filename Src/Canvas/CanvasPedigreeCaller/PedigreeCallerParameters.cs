using Newtonsoft.Json;
using Newtonsoft.Json.Converters;

namespace CanvasPedigreeCaller
{
    public enum CallerType
    {
        VariantCaller,
        HaplotypeVariantCaller
    }
    public class PedigreeCallerParameters
    {
        public int MaximumCopyNumber { get; set; } = 5;
        public int MaxAlleleNumber { get; set; } = 3;
        public int DefaultAlleleDensityThreshold { get; set; } = 1000;
        public double MaxQscore { get; set; } = 100.0;
        public int DefaultPerSegmentAlleleMaxCounts { get; set; } = 100;
        public int DefaultReadCountsThreshold { get; set; } = 4;
        public int MaxNumOffspringGenotypes { get; set; } = 500;
        public double DeNovoRate { get; set; } = 0.00001;
        public int MinimumCallSize { get; set; } = 2000;
        public int NumberOfTrimmedBins { get; set; } = 2;
        public int MaxCoreNumber { get; set; } = 30;
        [JsonProperty("DefaultCaller")]
        [JsonConverter(typeof(StringEnumConverter))]
        public CallerType DefaultCaller { get; set; } = CallerType.VariantCaller;
    }
}