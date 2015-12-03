using System.Collections.Generic;
using System.Linq;
using Illumina.Common;
using System;
using ProtoBuf;

namespace Isas.Shared
{
    // ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
    [Serializable]
    [ProtoContract]
    public class DescriptiveStats
    {
        [ProtoMember(1)]
        public float Max;
        [ProtoMember(2)]
        public float Mean;
        [ProtoMember(3)]
        public float Median;
        [ProtoMember(4)]
        public float Min;
        [ProtoMember(5)]
        public long NumberOfDataPoints;
        [ProtoMember(6)]
        public float Perc05;
        [ProtoMember(7)]
        public float Perc25;
        [ProtoMember(8)]
        public float Perc75;
        [ProtoMember(9)]
        public float Perc95;
        [ProtoMember(10)]
        public float Stdev;
        // ReSharper restore InconsistentNaming

        public DescriptiveStats()
        {
        }

        /// <summary>
        /// Generate descriptive stats from int array
        /// We make a copy first in case the caller assumes we don't modify the array
        /// </summary>
        public DescriptiveStats(int[] tempValues)
        {
            float[] values = Array.ConvertAll(tempValues, value => (float)value);
            GetStats(values, 1);
        }

        /// <summary>
        /// Generate descriptive stats from int array and optionally bin size
        /// We make a copy first in case the caller assumes we don't modify the array
        /// </summary>
        public DescriptiveStats(float[] tempValues, int binSize = 1)
        {
            float[] values = new float[tempValues.Length];
            Array.Copy(tempValues, values, tempValues.Length);
            GetStats(values, binSize);
        }

        private void GetStats(float[] values, int binSize)
        {
            if (values == null || values.Length == 0) return;
            Array.Sort(values);
            NumberOfDataPoints = values.Length * binSize;
            Mean = MathSupportFunctions.Mean(values) / binSize;
            Min = values.Min() / binSize;
            Max = values.Max() / binSize;
            Stdev = MathSupportFunctions.Stdev(values) / binSize;
            Median = MathSupportFunctions.PercentilePreSorted(values, 50) / binSize;
            Perc05 = MathSupportFunctions.PercentilePreSorted(values, 5) / binSize;
            Perc25 = MathSupportFunctions.PercentilePreSorted(values, 25) / binSize;
            Perc75 = MathSupportFunctions.PercentilePreSorted(values, 75) / binSize;
            Perc95 = MathSupportFunctions.PercentilePreSorted(values, 95) / binSize;
        }

        /// <summary>
        /// Generate descriptive stats from histogram
        /// </summary>
        public DescriptiveStats(IDictionary<int, long> hist)
        {
            //convert to the desired OrderedDictionary
            OrderedDictionary<long, long> orderedHist = new OrderedDictionary<long, long>();
            List<int> keys = hist.Keys.ToList();
            keys.Sort();
            foreach (int key in keys)
            {
                orderedHist.Add(key, hist[key]);
            }

            double[] percentilesIn = new double[0];
            long[] percentilesOut;
            InitializeFromHistogram(orderedHist, percentilesIn, out percentilesOut);
        }

        /// <summary>
        /// Generate descriptive stats from histogram
        /// </summary>
        public DescriptiveStats(OrderedDictionary<int, long> hist)
        {
            //convert to the desired OrderedDictionary
            OrderedDictionary<long, long> orderedHist = new OrderedDictionary<long, long>();
            List<int> keys = hist.Keys.ToList();
            keys.Sort();
            foreach (KeyValuePair<int, long> kvp in hist)
            {
                orderedHist.Add(kvp.Key, kvp.Value);
            }
            double[] percentilesIn = new double[0];
            long[] percentilesOut;
            InitializeFromHistogram(orderedHist, percentilesIn, out percentilesOut);
        }

        /// <summary>
        /// Generate descriptive stats from histogram
        /// </summary>
        public DescriptiveStats(OrderedDictionary<long, long> hist)
        {
            double[] percentilesIn = new double[0];
            long[] percentilesOut;
            InitializeFromHistogram(hist, percentilesIn, out percentilesOut);
        }

        /// <summary>
        /// Generate descriptive stats from histogram allowing for custom percentiles to be calcualted
        /// </summary>
        public DescriptiveStats(OrderedDictionary<long, long> hist, double[] extraPercentiles, out long[] calculatedPercentiles)
        {
            InitializeFromHistogram(hist, extraPercentiles, out calculatedPercentiles);
        }

        private void InitializeFromHistogram(OrderedDictionary<long, long> hist, double[] percentilesIn, out long[] percentilesOut)
        {
            percentilesOut = new long[percentilesIn.Length];
            Min = hist.Any() ? hist.First().Key : float.NaN;
            Max = hist.Any() ? hist.Last().Key : float.NaN;
            NumberOfDataPoints = hist.Values.Sum();

            long cumulativeCount = 0;
            double previousBinPercentile = 0;
            double currentBinPercentile = 0;

            //from http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Weighted_incremental_algorithm
            //D. H. D. West (1979). Communications of the ACM, 22, 9, 532-535: Updating Mean and Variance Estimates: An Improved Method
            long sumweight = 0;
            double mean = 0;
            double M2 = 0;
            foreach (KeyValuePair<long, long> bin in hist)
            {
                //handle mean and variance
                long x = bin.Key;
                long weight = bin.Value;
                long temp = weight + sumweight;
                double delta = x - mean;
                double R = delta * weight / temp;
                mean = mean + R;
                M2 = M2 + sumweight * delta * R;
                sumweight = temp;

                //handle percentiles
                cumulativeCount = sumweight;
                previousBinPercentile = currentBinPercentile;

                //percentile including values from this bin
                currentBinPercentile = ((double)cumulativeCount) / NumberOfDataPoints;

                //check if we have reached one of the percentiles
                if (0.05 > previousBinPercentile && 0.05 <= currentBinPercentile)
                {
                    Perc05 = bin.Key;
                }
                if (0.25 > previousBinPercentile && 0.25 <= currentBinPercentile)
                {
                    Perc25 = bin.Key;
                }
                if (0.5 > previousBinPercentile && 0.5 <= currentBinPercentile)
                {
                    Median = bin.Key;
                }
                if (0.75 > previousBinPercentile && 0.75 <= currentBinPercentile)
                {
                    Perc75 = bin.Key;
                }
                if (0.95 > previousBinPercentile && 0.95 <= currentBinPercentile)
                {
                    Perc95 = bin.Key;
                }

                //custom percentiles
                for (int percentileIndex = 0; percentileIndex < percentilesIn.Length; percentileIndex++)
                {
                    double percentile = percentilesIn[percentileIndex];
                    if (percentile > previousBinPercentile && percentile <= currentBinPercentile)
                    {
                        percentilesOut[percentileIndex] = bin.Key;
                    }
                }
            }
            Mean = (float)mean;
            if (NumberOfDataPoints > 1)
                Stdev = (float)Math.Sqrt(M2 / (NumberOfDataPoints - 1));
        }
    }
}