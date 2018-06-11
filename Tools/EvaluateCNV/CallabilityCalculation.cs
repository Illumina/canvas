using System.Collections.Generic;
using System.Linq;
using Illumina.Common;
using Isas.SequencingFiles;
using Isas.StatisticsProcessing.Processors.StatsProcessor;

namespace EvaluateCNV
{
    public class CallabilityCalculation : StatCalculation<CallabilityMetric, ReferenceInterval>
    {
        private readonly IReadOnlyDictionary<string, IEnumerable<Interval>> _intervals;

        public CallabilityCalculation(IReadOnlyDictionary<string, IEnumerable<Interval>> intervals)
        {
            _intervals = intervals;
        }

        public override CallabilityMetric Process(IEnumerable<ReferenceInterval> elements)
        {
            long totalBasePositions = _intervals.Sum(chromosomeKeyValuePair => chromosomeKeyValuePair.Value.Sum(interval => (long)interval.Length));
            long calledBasePositions = 0;
            foreach (var kvp in elements.GroupByAdjacent(bedEntry => bedEntry.Chromosome))
            {
                var chromosome = kvp.Key;
                var calledChromosome = kvp.Value.Select(element => element.Interval);
                var chromosomeIntervals = GetIntervalsForChromosome(chromosome);
                var overlap = new OverlapIntervals(chromosomeIntervals, calledChromosome);
                calledBasePositions += overlap.Sum(interval => interval.Length);
            }
            return new CallabilityMetric(calledBasePositions, totalBasePositions);
        }

        private IEnumerable<Interval> GetIntervalsForChromosome(string chromosome)
        {
            if (!_intervals.TryGetValue(chromosome, out var chromosomeIntervals)) chromosomeIntervals = Enumerable.Empty<Interval>();
            return chromosomeIntervals;
        }
    }
}