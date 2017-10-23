using System;
using System.Collections.Generic;
using System.Text;
using Xunit;
using static CanvasCommon.CanvasFilter;

namespace CanvasTest.CanvasCommon
{
    public class CanvasFilterTest
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
    }
}
