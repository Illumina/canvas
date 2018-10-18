using System.Collections.Generic;
using CanvasCommon;
using CanvasPartition.Models;

namespace CanvasPartition
{
    internal interface ISegmentationResultsProcessor
    {
        Dictionary<string, List<SegmentWithBins>> PostProcessSegments(GenomeSegmentationResults segmentationResults, PloidyInfo referencePloidy, Dictionary<string, List<SampleGenomicBin>> excludedIntervals, CoverageInfo coverageInfo);
    }
}