using System;
using System.Collections.Generic;
using CanvasCommon;
using Xunit;

namespace CanvasTest
{
    public class TestCanvasSomatic
    {
        [Fact]
        public void TestUsableSegments()
        {
            List<CanvasSegment> segments = new List<CanvasSegment>();
            int currentPosition = 1000;
            // Generate some segments.  Alternate between:
            // - Usable
            // - Too short
            // - Too few variants
            // - Too short + too few variants
            Random RNG = new Random();
            for (int index = 0; index < 100; index++)
            {
                int length = 100000;
                if (index % 2 == 1)
                {
                    length = 2000;
                }
                int variantCount = 999;
                if (index % 4 > 1) variantCount = 25;
                List<SampleGenomicBin> counts = new List<SampleGenomicBin>();
                for (int countIndex = 0; countIndex < length / 100; countIndex++) counts.Add(new SampleGenomicBin("1", 2, 3, RNG.Next(1000)));
                CanvasSegment segment = new CanvasSegment("chr1", currentPosition, currentPosition + length, counts);
                for (int varIndex = 0; varIndex < variantCount; varIndex++)
                {
                    segment.Balleles.Add(new Ballele(0, RNG.Next(), 100, 50, 50));
                }
                segments.Add(segment);
            }
            var usable = CanvasSomaticCaller.SomaticCaller.GetUsableSegmentsForModeling(segments, false, 50, null);
            Assert.Equal(50, usable.Count);
        }
    }
}
