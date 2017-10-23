using System.Collections.Generic;
using System.Linq;

namespace CanvasCommon
{
    public static class CanvasFilter
    {
        public const int SegmantSizeCutoff = 10000;
        public const string Pass = "PASS";
        public static List<string> PassedFilter = new List<string> {Pass};
        public const string AllSampleFiltersFailed = "FailedFT";
        public const string SampleFilterFailed = "FT";
        public static string FormatCnvSizeWithSuffix(int size)
        {
            if (size >= 1000000)
                return FormatCnvSizeWithSuffix(size / 1000000) + "MB";
            if (size >= 1000)
                return FormatCnvSizeWithSuffix(size / 1000) + "KB";
            return size.ToString();
        }

        // The record level fiter keywords in Canvas are always derived from sample level filters
        public static IEnumerable<string> GetRecordLevelFilter(IReadOnlyCollection<List<string>> sampleLevelFiters)
        {
            var keywords = sampleLevelFiters.SelectMany(x => x).GroupBy(x => x);
            var failedRecordLevelFilter = new List<string> { AllSampleFiltersFailed };
            var nSamples = sampleLevelFiters.Count;
            foreach (var keyword in keywords)
            {
                if (keyword.Key == Pass) return PassedFilter; // At least one sample filter passed
                if (keyword.Count() == nSamples) failedRecordLevelFilter.Add(keyword.Key); // Failed for all samples
            }
            return failedRecordLevelFilter;
        }

        public static string ToString(IEnumerable<string> filter)
        {
            return string.Join(";", filter);
        }
    }
}