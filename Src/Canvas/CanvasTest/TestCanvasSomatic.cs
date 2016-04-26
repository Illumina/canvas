using System;
using System.Collections.Generic;
using CanvasCommon;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace CanvasTest
{
    [TestClass]
    public class TestCanvasSomatic
    {
        [TestMethod]
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
                List<float> counts = new List<float>();
                for (int countIndex = 0; countIndex < length / 100; countIndex++) counts.Add(RNG.Next(1000));
                CanvasSegment segment = new CanvasSegment("chr1", currentPosition, currentPosition + length, counts);
                for (int varIndex = 0; varIndex < variantCount; varIndex++)
                {
                    segment.VariantFrequencies.Add(RNG.Next());
                }
                segments.Add(segment);
            }
            var usable = CanvasSomaticCaller.SomaticCaller.GetUsableSegmentsForModeling(segments, false, 50);
            Assert.AreEqual(50, usable.Count);
        }
    }
}
