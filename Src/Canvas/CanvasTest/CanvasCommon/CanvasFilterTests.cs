using System.Collections.Generic;
using CanvasCommon;
using Xunit;

namespace CanvasTest
{
    public class CanvasFilterTests
    {
        [Fact]
        public void FormatCnvSizeWithSuffix_10kb_conversion()
        {
            const int size = 10000;
            var formattedSize = CanvasFilter.FormatCnvSizeWithSuffix(size);
            Assert.Equal("10KB", formattedSize);
        }

        [Fact]
        public void FormatCnvSizeWithSuffix_6mb_conversion()
        {
            const int size = 6000000;
            var formattedSize = CanvasFilter.FormatCnvSizeWithSuffix(size);
            Assert.Equal("6MB", formattedSize);
        }

        [Fact]
        public void FormatCnvSizeWithSuffix_500bp_conversion()
        {
            const int size = 500;
            var formattedSize = CanvasFilter.FormatCnvSizeWithSuffix(size);
            Assert.Equal("500", formattedSize);
        }

        [Fact]
        public void ToVcfString_test()
        {
            var canvasFilter = CanvasFilter.Create(new[] { "First", "Second", "Third" });
            Assert.Equal("First;Second;Third", canvasFilter.ToVcfString());
        }

        [Fact]
        public void Create_Filter_w_pass_tag_ignored()
        {
            var filterTags = new[] { "First", "Second", CanvasFilter.Pass };
            var canvasFilter = CanvasFilter.Create(filterTags);
            Assert.Equal("First;Second" , canvasFilter.ToVcfString());
        }

        [Fact]
        public void GetRecordLevelFilterFromSampleFiltersOnly_single_sample_pass()
        {
            var sampleLevelFilters = new[] { CanvasFilter.PassFilter };
            var recordLevelFilter = CanvasFilter.GetRecordLevelFilterFromSampleFiltersOnly(sampleLevelFilters);
            Assert.Equal("PASS", recordLevelFilter.ToVcfString());
        }

        [Fact]
        public void GetRecordLevelFilterFromSampleFiltersOnly_multi_sample_all_pass()
        {
            var sampleLevelFilters = new[] { CanvasFilter.PassFilter, CanvasFilter.PassFilter, CanvasFilter.PassFilter };
            var recordLevelFilter = CanvasFilter.GetRecordLevelFilterFromSampleFiltersOnly(sampleLevelFilters);
            Assert.Equal("PASS", recordLevelFilter.ToVcfString());
        }

        [Fact]
        public void GetRecordLevelFilterFromSampleFiltersOnly_multi_sample_one_pass()
        {
            var sampleLevelFilters = new[]
            {
                CanvasFilter.PassFilter,
                CanvasFilter.Create(new [] { "Failed1" }),
                CanvasFilter.Create(new [] { "Failed2" }),
            };
            var recordLevelFilter = CanvasFilter.GetRecordLevelFilterFromSampleFiltersOnly(sampleLevelFilters);
            Assert.Equal("PASS", recordLevelFilter.ToVcfString());
        }

        [Fact]
        public void GetRecordLevelFilter_multi_sample_failed_different_filters()
        {
            var sampleLevelFilters = new[]
            {
                CanvasFilter.Create(new [] { "Failed1" }),
                CanvasFilter.Create(new [] { "Failed2" }),
                CanvasFilter.Create(new [] { "Failed3" })
            };
            var recordLevelFilter = CanvasFilter.GetRecordLevelFilterFromSampleFiltersOnly(sampleLevelFilters);
            Assert.Equal("FailedFT", recordLevelFilter.ToVcfString());
        }

        [Fact]
        public void GetRecordLevelFilter_multi_sample_failed_common_filter()
        {
            var sampleLevelFilters = new[]
            {
                CanvasFilter.Create(new [] { "Failed1" }),
                CanvasFilter.Create(new [] { "Failed1", "Failed2" }),
                CanvasFilter.Create(new [] { "Failed1", "Failed3" })
            };
            var recordLevelFilter = CanvasFilter.GetRecordLevelFilterFromSampleFiltersOnly(sampleLevelFilters);
            Assert.Equal("FailedFT", recordLevelFilter.ToVcfString());
        }

        [Fact]
        public void UpdateRecordLevelFilter_w_MinQual()
        {
            var sampleLevelFilters = new[]
            {
                CanvasFilter.Create(new [] { "Failed1" }),
                CanvasFilter.Create(new [] { "Failed1", "Failed2" }),
                CanvasFilter.Create(new [] { "Failed1", "Failed3" })
            };
            var recordLevelFitler = CanvasFilter.Create(new[] { "MinQual" });
            var updatedrecordLevelFilter = CanvasFilter.UpdateRecordLevelFilter(recordLevelFitler, sampleLevelFilters);
            Assert.Equal("MinQual;FailedFT", updatedrecordLevelFilter.ToVcfString());
        }
    }
}
