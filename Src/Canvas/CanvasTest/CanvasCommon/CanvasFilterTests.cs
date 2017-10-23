using System;
using System.Collections.Generic;
using System.Text;
using CanvasCommon;
using Xunit;
using static CanvasCommon.CanvasFilter;

namespace CanvasTest.CanvasCommon
{
    public class CanvasFilterTests
    {
        [Fact]
        public void FormatCnvSizeWithSuffix_10kb_conversion()
        {
            var size = 10000;
            var formattedSize = FormatCnvSizeWithSuffix(size);
            Assert.Equal("10KB", formattedSize);
        }

        [Fact]
        public void FormatCnvSizeWithSuffix_6mb_conversion()
        {
            var size = 6000000;
            var formattedSize = FormatCnvSizeWithSuffix(size);
            Assert.Equal("6MB", formattedSize);
        }

        [Fact]
        public void FormatCnvSizeWithSuffix_500bp_conversion()
        {
            var size = 500;
            var formattedSize = FormatCnvSizeWithSuffix(size);
            Assert.Equal("500", formattedSize);
        }

        [Fact]
        public void ToString_test()
        {
            var recordLevelFilter = new List<string> { "First", "Second", "Third" };
            Assert.Equal("First;Second;Third", CanvasFilter.ToString(recordLevelFilter));
        }

        [Fact]
        public void GetRecordLevelFilter_single_sample_pass()
        {
            var sampleLevelFilters = new List<List<string>> { PassedFilter };
            var recordLevelFilter = GetRecordLevelFilter(sampleLevelFilters);
            Assert.Equal(new List<string> { "PASS" }, recordLevelFilter);
        }

        [Fact]
        public void GetRecordLevelFilter_multi_sample_all_pass()
        {
            var sampleLevelFilters = new List<List<string>> { PassedFilter, PassedFilter, PassedFilter };
            var recordLevelFilter = GetRecordLevelFilter(sampleLevelFilters);
            Assert.Equal(new List<string> { "PASS" }, recordLevelFilter);
        }

        [Fact]
        public void GetRecordLevelFilter_multi_sample_one_pass()
        {
            var sampleLevelFilters = new List<List<string>>();
            sampleLevelFilters.Add(PassedFilter);
            sampleLevelFilters.Add(new List<string> { "Failed1" });
            sampleLevelFilters.Add(new List<string> { "Failed2" });
            var recordLevelFilter = GetRecordLevelFilter(sampleLevelFilters);
            Assert.Equal(new List<string> { "PASS" }, recordLevelFilter);
        }

        [Fact]
        public void GetRecordLevelFilter_multi_sample_failed_different_filters()
        {
            var sampleLevelFilters = new List<List<string>>();
            sampleLevelFilters.Add(new List<string> { "Failed1" });
            sampleLevelFilters.Add(new List<string> { "Failed2" });
            sampleLevelFilters.Add(new List<string> { "Failed3" });
            var recordLevelFilter = GetRecordLevelFilter(sampleLevelFilters);
            Assert.Equal(new List<string> { "FailedFT" }, recordLevelFilter);
        }

        [Fact]
        public void GetRecordLevelFilter_multi_sample_failed_common_filter()
        {
            var sampleLevelFilters = new List<List<string>>();
            sampleLevelFilters.Add(new List<string> { "Failed1" });
            sampleLevelFilters.Add(new List<string> { "Failed1", "Failed2" });
            sampleLevelFilters.Add(new List<string> { "Failed1", "Failed3" });
            var recordLevelFilter = GetRecordLevelFilter(sampleLevelFilters);
            Assert.Equal(new List<string> { "FailedFT", "Failed1" }, recordLevelFilter);
        }
    }
}
