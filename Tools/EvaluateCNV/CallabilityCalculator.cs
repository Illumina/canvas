using System.Collections.Generic;
using Isas.SequencingFiles;
using Isas.StatisticsProcessing.Processors.StatsProcessor;
using Isas.StatisticsProcessing.SummaryMetrics;

namespace EvaluateCNV
{
    public class CallabilityCalculator : StatCalculator<CallabilityMetric, ReferenceInterval>
    {
        private readonly IReadOnlyDictionary<string, IEnumerable<Interval>> _intervals;

        public CallabilityCalculator(IReadOnlyDictionary<string, IEnumerable<Interval>> intervals)
        {
            _intervals = intervals;
        }

        public override StatCalculation<CallabilityMetric, ReferenceInterval> NewCalculation()
        {
            return new CallabilityCalculation(_intervals);
        }
    }
}