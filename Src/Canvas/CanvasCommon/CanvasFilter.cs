using System;
using System.Collections.Generic;
using System.Linq;
using Illumina.Common;

namespace CanvasCommon
{
    public class CanvasFilter
    {

        public const int SegmantSizeCutoff = 10000;
        public const string Pass = "PASS";
        public static List<string> PassedFilter = new List<string> { Pass };
        public const string AllSampleFiltersFailed = "FailedFT";
        public const string SampleFilterFailed = "FT";
        public IReadOnlyList<string> SetOfFailedFilters { get; }

        private CanvasFilter(IReadOnlyList<string> setOfFailedFilters)
        {
            SetOfFailedFilters = setOfFailedFilters;
        }

        public static CanvasFilter PassFilter = Create(Enumerable.Empty<string>());

        public static CanvasFilter Create(IEnumerable<string> setOfFailedFilters)
        {
            var allFailedFilters = new List<string>();
            foreach (var failedFilters in setOfFailedFilters)
            {
                if (failedFilters == Pass)
                    throw new IlluminaException($"{nameof(failedFilters)} cannot contain a filter named '{Pass}'");

                allFailedFilters.Add(failedFilters);
            }
            return new CanvasFilter(allFailedFilters.ToReadOnlyList());
        }

        public string ToVcfString() => IsPass ? Pass : string.Join(";", SetOfFailedFilters);

        public bool IsPass => SetOfFailedFilters.Count == 0;

        public static CanvasFilter GetSharedFilters(IReadOnlyList<CanvasFilter> sampleFilters)
        {
            var sharedFilters = sampleFilters.First().SetOfFailedFilters;
            foreach (var sampleFilter in sampleFilters.Skip(1))
            {
                if (sharedFilters.Empty()) return PassFilter;
                sharedFilters = sharedFilters.Intersect(sampleFilter.SetOfFailedFilters).ToList();
            }
            return Create(sharedFilters);
        }

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

        public static string ToString(IEnumerable<string> filter) => string.Join(";", filter);
    }
}
