
using System.Collections.Generic;
using System.Linq;
using Illumina.Common;

namespace CanvasCommon
{
    public class CanvasFilter
    {

        public const int SegmentSizeCutoff = 10000;
        public const string Pass = "PASS";
        public const string AllSampleFailedTag = "FailedFT";
        public const string SampleFilterFailed = "FT";
        public static readonly CanvasFilter PassFilter = Create(Enumerable.Empty<string>());
        private IReadOnlyList<string> FailedFilterTags { get; }

        private CanvasFilter(IReadOnlyList<string> failedFilterTags)
        {
            FailedFilterTags = failedFilterTags;
        }

        public static CanvasFilter Create(IEnumerable<string> filterTags)
        {
            var allFailedFilterTags = filterTags.Where(failedFilterTag => failedFilterTag != Pass).ToList();
            return new CanvasFilter(allFailedFilterTags.ToReadOnlyList());
        }

        public string ToVcfString() => IsPass ? Pass : string.Join(";", FailedFilterTags);

        public bool IsPass => FailedFilterTags.Count == 0;

        private static CanvasFilter GetSharedFilter(IReadOnlyList<CanvasFilter> sampleFilters)
        {
            if (sampleFilters.First().IsPass) return PassFilter;
            var sharedFilterTags = sampleFilters.First().FailedFilterTags.ToList();
            foreach (var sampleFilter in sampleFilters.Skip(1))
            {
                if (sampleFilter.IsPass) return PassFilter;
                sharedFilterTags = sharedFilterTags.Intersect(sampleFilter.FailedFilterTags).ToList();
            }
            sharedFilterTags.Add(AllSampleFailedTag);
            return Create(sharedFilterTags);
        }

        public static CanvasFilter UpdateRecordLevelFilter(CanvasFilter recordLevelFilter, IReadOnlyList<CanvasFilter> sampleFilters) => Create(recordLevelFilter.FailedFilterTags.Concat(GetSharedFilter(sampleFilters).FailedFilterTags));

        public static CanvasFilter GetRecordLevelFilterFromSampleFiltersOnly(
            IReadOnlyList<CanvasFilter> sampleFilters) => UpdateRecordLevelFilter(PassFilter, sampleFilters);

        public static string FormatCnvSizeWithSuffix(int size)
        {
            if (size >= 1000000)
                return FormatCnvSizeWithSuffix(size / 1000000) + "MB";
            if (size >= 1000)
                return FormatCnvSizeWithSuffix(size / 1000) + "KB";
            return size.ToString();
        }
    }
}
