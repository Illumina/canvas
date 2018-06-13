using System.Collections.Generic;
using Isas.StatisticsProcessing.SummaryMetrics;

namespace EvaluateCNV
{
    public class CallabilityMetric : ISummaryMetricProvider
    {
        public long CalledBases { get; }
        public long TotalBases { get; }
        public double Callability => ((double)CalledBases) / TotalBases;

        public CallabilityMetric(long calledBases, long totalBases)
        {
            CalledBases = calledBases;
            TotalBases = totalBases;
        }
        
        public IEnumerable<SummaryMetricBase> GetMetrics()
        {
            var totalBasePositionsMetric = new SummaryMetric<long>("Total base positions")
            {
                Value = TotalBases
            };
            var calledBasePositionsMetric = new SummaryMetric<long>("Called base positions")
            {
                Value = CalledBases
            };
            var callabilityPercent = new SummaryMetricPercent<long, long>(calledBasePositionsMetric, totalBasePositionsMetric, "Percent callability");
            yield return totalBasePositionsMetric;
            yield return calledBasePositionsMetric;
            yield return callabilityPercent;
        }
    }
}