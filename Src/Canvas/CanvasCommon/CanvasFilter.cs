using System.Collections.Generic;

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
            if (size >= 100000000)
                return FormatCnvSizeWithSuffix(size / 1000000) + "MB";
            if (size >= 100000)
                return FormatCnvSizeWithSuffix(size / 1000) + "KB";
            if (size >= 10000)
            {
                return (size / 1000D).ToString("0.#") + "KB";
            }
            return size.ToString();
        }
    }
}