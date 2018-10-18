using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Illumina.Common.Collections;

namespace CanvasPartition
{
    internal class GenomeSegmentationResults
    {
        public IDictionary<string, SegmentationInput.Segment[]> SegmentByChr;

        public GenomeSegmentationResults(IDictionary<string, SegmentationInput.Segment[]> segmentByChr)
        {
            this.SegmentByChr = segmentByChr;
        }

        public static GenomeSegmentationResults SplitOverlappingSegments(
            List<GenomeSegmentationResults> sampleSegmentationResults)
        {
            if (sampleSegmentationResults.Count == 1) return sampleSegmentationResults.Single();

            var result = new ConcurrentDictionary<string, SegmentationInput.Segment[]>();
            Parallel.ForEach(sampleSegmentationResults.First().SegmentByChr.Keys, chr =>
            {
                result[chr] = SplitOverlappingSegments(sampleSegmentationResults
                    .Select(segmentation => segmentation.SegmentByChr[chr]).ToList()).ToArray();
            });
            return new GenomeSegmentationResults(result);
        }

        private static IEnumerable<SegmentationInput.Segment> SplitOverlappingSegments(List<SegmentationInput.Segment[]> sampleSegments)
        {
            var starts = MergeEnumerator.Merge(sampleSegments.Select(segments => segments.Select(segment => segment.start)));
            var ends = MergeEnumerator.Merge(sampleSegments.Select(segments => segments.Select(segment => segment.end)));

            var partitions = MergeEnumerator.Merge(new[]
            {
                starts.Select(start => (Position: start, IsStart: true)),
                ends.Select(end => (Position: end, IsStart: false))
            }, (position1, position2) => position1.CompareTo(position2));

            var numberOverlappingSegments = 0;
            uint currentPosition = 0;
            foreach (var (position, isStart) in partitions)
            {
                if (numberOverlappingSegments > 0 && currentPosition != position)
                {
                    yield return new SegmentationInput.Segment { start = currentPosition, end = position };
                }

                currentPosition = position;
                numberOverlappingSegments += isStart ? 1 : -1;
            }
        }
    }
}