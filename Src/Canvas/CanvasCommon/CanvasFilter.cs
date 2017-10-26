using System;
using System.Collections.Generic;
using System.Linq;
using Illumina.Common;

namespace CanvasCommon
{
    public class CanvasFilter
    {
        public IReadOnlyList<string> FailedFilters { get; }

        private CanvasFilter(IReadOnlyList<string> failedFilters)
        {
            FailedFilters = failedFilters;
        }

        public static CanvasFilter PassFilter = Create(Enumerable.Empty<string>());

        public static CanvasFilter Create(IEnumerable<string> failedFilters)
        {
            var filters = new List<string>();
            foreach (var filter in failedFilters)
            {
                if (filter == Pass)
                    throw new IlluminaException($"{nameof(failedFilters)} cannot contain a filter named '{Pass}'");
                if (filters.Contains(filter))
                    throw new ArgumentException($"{nameof(failedFilters)} contains duplicate filter {filter}");
                filters.Add(filter);
            }
            return new CanvasFilter(filters.ToReadOnlyList());
        }

        public string ToVcfString()
        {
            if (IsPass)
                return Pass;
            return string.Join(";", FailedFilters);
        }

        public bool IsPass => FailedFilters.Count == 0;

        public static CanvasFilter GetSharedFilters(IReadOnlyList<CanvasFilter> sampleFilters)
        {
            var sharedFilters = sampleFilters.First().FailedFilters;
            foreach (var sampleFilter in sampleFilters.Skip(1))
            {
                if (sharedFilters.Empty()) return PassFilter;
                sharedFilters = sharedFilters.Intersect(sampleFilter.FailedFilters).ToList();
            }
            return Create(sharedFilters);
        }

        public const int SegmantSizeCutoff = 10000;
        public const string Pass = "PASS";
        public static List<string> PassedFilter = new List<string> { Pass };
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

        public static string ToString(IEnumerable<string> filter) => string.Join(";", filter);
    }
}
