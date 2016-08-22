using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;

namespace CanvasCommon
{

    public enum CanvasCoverageMode
    {
        Binary = 0, // Count 0 or 1 hits per kmer; not used any more!
        TruncatedDynamicRange, // was 3
        GCContentWeighted, // was 5
        Fragment // Count fragments not unique 35-mers
    }

    public enum CanvasNormalizeMode
    {
        WeightedAverage, // Weighted average of all control samples
        BestLR2, // Use sum of log ratio^2 to find the best control sample
        PCA // Project a sample to the PCA axes to obtain the control sample
    }

    public enum CanvasGCNormalizationMode
    {
        MedianByGC, // Smooth by (global median) * count / (median at %GC)
        LOESS // Smooth by (global median) + count - (fitted count). Counts are transformed by sqrt.
    }

    public enum CanvasSomaticClusteringMode
    {
        GaussianMixture, 
        Density 
    }

    public static class Utilities
    {
        #region Members
        static private Dictionary<int, double> CachedLogFactorial;
        #endregion


        /// <summary>
        /// utility method for transforming the filename extension (e.g. from .genome.vcf.gz to .vcf)
        /// Using String.Replace method will replace all occurences while this method only replaces the occurence at the end
        /// </summary>
        public static string ReplaceFileNameExtension(this string path, string oldSuffix, string newSuffix)
        {
            if (!path.EndsWith(oldSuffix)) return path;
            return path.Substring(0, path.Length - oldSuffix.Length) + newSuffix;
        }

        static public string GetCoverageAndVariantFrequencyOutputPath(string outputVcfPath)
        {
            string coveragePath = outputVcfPath;
            if (outputVcfPath.EndsWith(".vcf.gz"))
                coveragePath = coveragePath.ReplaceFileNameExtension(".vcf.gz", "");
            else
                coveragePath.ReplaceFileNameExtension(".vcf", "");
            return coveragePath + ".CoverageAndVariantFrequency.txt";
        }


        static public CanvasCoverageMode ParseCanvasCoverageMode(string mode)
        {
            switch (mode.ToLowerInvariant().Trim())
            {
                case "0":
                case "binary":
                    return CanvasCoverageMode.Binary;
                case "truncateddynamicrange":
                case "3":
                    return CanvasCoverageMode.TruncatedDynamicRange;
                case "5":
                case "gccontentweighted":
                    return CanvasCoverageMode.GCContentWeighted;
                case "fragment":
                    return CanvasCoverageMode.Fragment;
                default:
                    throw new Exception(string.Format("Invalid canvas coverage mode '{0}'", mode));
            }
        }

        static public CanvasNormalizeMode ParseCanvasNormalizeMode(string mode)
        {
            switch (mode.ToLowerInvariant().Trim())
            {
                case "weightedaverage":
                    return CanvasNormalizeMode.WeightedAverage;
                case "bestlr2":
                    return CanvasNormalizeMode.BestLR2;
                case "pca":
                    return CanvasNormalizeMode.PCA;
                default:
                    throw new Exception(string.Format("Invalid CanvasNormalize mode '{0}'", mode));
            }
        }

        static public CanvasGCNormalizationMode ParseCanvasGCNormalizationMode(string mode)
        {
            switch (mode.ToLowerInvariant().Trim())
            {
                case "medianbygc":
                    return CanvasGCNormalizationMode.MedianByGC;
                case "loess":
                    return CanvasGCNormalizationMode.LOESS;
                default:
                    throw new Exception(string.Format("Invalid CanvasClean mode '{0}'", mode));
            }
        }

        static public CanvasSomaticClusteringMode ParseCanvasSomaticClusteringMode(string mode)
        {
            switch (mode.ToLowerInvariant().Trim())
            {
                case "gaussianmixture":
                    return CanvasSomaticClusteringMode.GaussianMixture;
                case "density":
                    return CanvasSomaticClusteringMode.Density;
                default:
                    throw new Exception(string.Format("Invalid CanvasSomatic clustering mode '{0}'", mode));
            }
        }

        static public void LogCommandLine(string[] arguments)
        {
            Console.WriteLine(">>>Command-line arguments:");
            foreach (string arg in arguments) Console.Write("{0} ", arg);
            Console.WriteLine();
        }

        static public double Mean(int[] x)
        {
            double sum = 0;

            for (int i = 0; i < x.Length; i++)
                sum += x[i];

            return sum / x.Length;
        }

        // calculate mean of Int16 values greaer than zero
        static public Int16 NonZeroMean(Int16[] x)
        {
            long sum = 0;
            long counter = 0;

            for (int i = 0; i < x.Length; i++)
            {
                if (x[i] > 0)
                {
                    sum += x[i];
                    counter++;
                }
            }
            if (counter == 0) return 0;
            return Convert.ToInt16(sum / counter);
        }

        // calculate weighted mean
        static public double WeightedMean(List<double> x, List<double> weight)
        {
            double sum = 0;
            double counter = 0;

            for (int i = 0; i < x.Count; i++)
            {
                sum += x[i] * weight[i];
                counter += weight[i];
            }

            if (counter == 0) return 0;
            else return (sum / counter);
        }

        // calculate mean of non-zero byte values
        static public Int16 NonZeroMean(byte[] x)
        {
            long sum = 0;
            long counter = 0;

            for (int i = 0; i < x.Length; i++)
            {
                if (x[i] > 0)
                {
                    sum += x[i];
                    counter++;
                }
            }
            if (counter == 0) return 0;
            return Convert.ToInt16(sum / counter);
        }


        static public double Mean(float[] x)
        {
            double sum = 0;

            for (int i = 0; i < x.Length; i++)
                sum += x[i];

            return sum / x.Length;
        }

        // estimate Mean of an arbitrary interval 
        static public double Mean(double[] x, int start = 0, int end = 0)
        {
            if (start == 0 && end == 0)
            {
                end = x.Length;
            }

            double sum = 0;

            for (int i = start; i < end; i++)
                sum += x[i];

            return sum / (end - start);
        }

        static public double StandardDeviation(int[] x)
        {
            double mu = Mean(x);
            double sum = 0;

            for (int i = 0; i < x.Length; i++)
            {
                double diff = x[i] - mu;
                sum += diff * diff;
            }
            return Math.Sqrt(sum / (x.Length - 1));
        }

        static public double StandardDeviation(float[] x)
        {
            double mu = Mean(x);
            double sum = 0;

            for (int i = 0; i < x.Length; i++)
            {
                double diff = x[i] - mu;
                sum += diff * diff;
            }
            return Math.Sqrt(sum / (x.Length - 1));
        }

        // estimate Standard Deviation of an arbitrary interval 
        static public double StandardDeviation(double[] x, int start = 0, int end = 0)
        {
            if (start == 0 && end == 0)
            {
                end = x.Length;
            }

            double mu = Mean(x, start, end);
            double sum = 0;

            for (int i = start; i < end; i++)
            {
                double diff = x[i] - mu;
                sum += diff * diff;
            }
            return Math.Sqrt(sum / (end - start - 1));
        }

        // estimate Standard Deviation of double
        static public double StandardDeviation(double[] x)
        {

            double mu = Mean(x);
            double sum = 0;

            for (int i = 0; i < x.Length; i++)
            {
                double diff = x[i] - mu;
                sum += diff * diff;
            }
            return Math.Sqrt(sum / (x.Length - 1));
        }

        // estimate Standard Deviation of double list
        static public double StandardDeviation(List<double> x)
        {

            double mu = x.Average();
            double sum = 0;

            for (int i = 0; i < x.Count; i++)
            {
                double diff = x[i] - mu;
                sum += diff * diff;
            }
            return Math.Sqrt(sum / (x.Count - 1));
        }

        /// <summary>
        /// Calculates the median of a list.  Does not re-order the original list.
        /// </summary>
        /// <param name="x">List of doubles.</param>
        /// <returns>Median of x.</returns>
        public static double Median(IEnumerable<double> x)
        {

            List<double> sorted = new List<double>(x.Count());

            foreach (double value in x)
            {
                sorted.Add(value);
            }

            sorted.Sort();

            int n = sorted.Count;

            double median = 0;

            if (n % 2 == 0)
            {
                double val1 = sorted[(n / 2) - 1];
                double val2 = sorted[(n / 2)];
                median = (val1 + val2) / 2;
            }
            else
            {
                median = sorted[(n / 2)];
            }

            return median;
        }


        /// <summary>
        /// Calculates first and third quartiles of the input data vector
        /// This method uses the following logic to process odd and even vector lengths (to make robust against small re-sequencing panels):
        /// If there are an even number of data points:
        ///    - Use the median to divide the ordered data set into two halves. 
        ///    - The lower quartile value is the median of the lower half of the data. 
        ///    - The upper quartile value is the median of the upper half of the data.
        /// If there are (4n+1) data points:
        ///    - The lower quartile is 25% of the nth data value plus 75% of the (n+1)th data value.
        ///    - The upper quartile is 75% of the (3n+1)th data point plus 25% of the (3n+2)th data point.
        ///If there are (4n+3) data points:
        ///    - The lower quartile is 75% of the (n+1)th data value plus 25% of the (n+2)th data value.
        ///    - The upper quartile is 25% of the (3n+2)th data point plus 75% of the (3n+3)th data point.
        /// <param name="x">List of doubles.</param>
        /// <returns>First, second and third quartiles of x.</returns>
        static public Tuple<float, float, float> Quartiles(List<float> x)
        {
            List<float> sorted = new List<float>(x.Count);

            for (int i = 0; i < x.Count; i++)
            {
                sorted.Add(x[i]);
            }

            sorted.Sort();

            int iSize = sorted.Count;
            int iMid = iSize / 2; //this is the mid from a zero based index, eg mid of 7 = 3;

            float fQ1 = 0;
            float fQ2 = 0;
            float fQ3 = 0;

            // even number of points
            if (iSize % 2 == 0)
            {
                fQ2 = (sorted[iMid - 1] + sorted[iMid]) / 2;
                int iMidMid = iMid / 2;
                //easy split 
                if (iMid % 2 == 0)
                {
                    fQ1 = (sorted[iMidMid - 1] + sorted[iMidMid]) / 2;
                    fQ3 = (sorted[iMid + iMidMid - 1] + sorted[iMid + iMidMid]) / 2;
                }
                else
                {
                    fQ1 = sorted[iMidMid];
                    fQ3 = sorted[iMidMid + iMid];
                }
            }

            // odd number so the median is just the midpoint in the array
            else
            {
                fQ2 = sorted[iMid];
                // (4n-1) points
                if ((iSize - 1) % 4 == 0)
                {
                    int n = (iSize - 1) / 4;
                    fQ1 = (sorted[n - 1] * 0.25f) + (sorted[n] * 0.75f);
                    fQ3 = (sorted[3 * n] * 0.75f) + (sorted[3 * n + 1] * 0.25f);
                }
                // (4n-3) points   
                else if ((iSize - 3) % 4 == 0)
                {
                    int n = (iSize - 3) / 4;

                    fQ1 = (sorted[n] * 0.75f) + (sorted[n + 1] * 0.25f);
                    fQ3 = (sorted[3 * n + 1] * 0.25f) + (sorted[3 * n + 2] * 0.75f);
                }
            }

            return new Tuple<float, float, float>(fQ1, fQ2, fQ3);
        }

        /// <summary>
        /// Calculates the median of an arbitrary interval in a  list.  Does not re-order the original list.
        /// </summary>
        /// <param name="x">List of doubles.</param>
        /// <param name="start">First element to calculate the median.</param>
        /// <param name="end">(Exclusive) Last element to calculate the median.</param>
        /// <returns>Median of x.</returns>
        public static double Median(IEnumerable<double> x, int start = 0, int? end = null)
        {
            if (start < 0)
                throw new ArgumentException("start cannot be negative");
            if (end.HasValue && end.Value <= start)
                throw new ArgumentException("end must be greater than start");

            double[] sorted;
            if (!end.HasValue)
                sorted = x.Skip(start).ToArray();
            else
                sorted = x.Skip(start).Take(end.Value - start).ToArray();
            Array.Sort(sorted);

            int n = sorted.Length;
            double median = 0;

            if (n % 2 == 0)
            {
                double val1 = sorted[(n / 2) - 1];
                double val2 = sorted[(n / 2)];
                median = (val1 + val2) / 2;
            }
            else
            {
                median = sorted[(n / 2)];
            }

            return median;
        }

        /// <summary>
        /// Calculates the median absolute deviation (MAD) of a list.  Does not re-order the original list.
        /// </summary>
        /// <param name="x">List of doubles.</param>
        /// <param name="start">First element to calculate the MAD.</param>
        /// <param name="end">(Exclusive) Last element to calculate the MAD.</param>
        /// <returns>Median absolute deviation of x.</returns>
        public static double Mad(IEnumerable<double> x, int start = 0, int? end = null)
        {
            if (start < 0)
                throw new ArgumentException("start cannot be negative");
            if (end.HasValue && end.Value <= start)
                throw new ArgumentException("end must be greater than start");

            double median = Median(x, start, end);
            var xSlice = !end.HasValue ? x.Skip(start) : x.Skip(start).Take(end.Value - start);
            double[] diffs = xSlice.Select(xi => Math.Abs(xi - median)).ToArray();

            return Median(diffs);
        }

        /// <summary>
        /// Calculates the median of a list.  Does not re-order the original list.
        /// </summary>
        /// <param name="x">List of doubles.</param>
        /// <returns>Median of x.</returns>
        public static double Median(List<float> x)
        {
            int n = x.Count;

            List<float> sorted = new List<float>(n);

            foreach (float v in x)
                sorted.Add(v);

            sorted.Sort();

            double median = 0;

            if (n % 2 == 0)
            {
                double val1 = sorted[(n / 2) - 1];
                double val2 = sorted[(n / 2)];
                median = (val1 + val2) / 2;
            }
            else
            {
                median = sorted[(n / 2)];
            }

            return median;
        }

        /// <summary>
        /// Computes the median value of a list of ints.
        /// </summary>
        /// <param name="x">List of ints to summarize.</param>
        /// <returns></returns>
        public static int Median(IEnumerable<int> x)
        {
            int n = x.Count();
            if (n == 0) return 0;
            List<int> sorted = new List<int>(n);

            foreach (int v in x)
                sorted.Add(v);

            sorted.Sort();

            int median = 0;

            if (n % 2 == 0)
            {
                int val1 = sorted[(n / 2) - 1];
                int val2 = sorted[(n / 2)];
                median = (val1 + val2) / 2;
            }
            else
            {
                median = sorted[(n / 2)];
            }

            return median;

        }

        /// <summary>
        /// Calculate the weighted quantiles. The weights must be non-negative
        /// </summary>
        /// <param name="x">List of (number, weight) tuples</param>
        /// <param name="probs"></param>
        /// <returns></returns>
        public static double[] WeightedQuantiles(List<Tuple<float, float>> x, List<float> probs) 
        {
            if (x.Any(t => t.Item2 < 0)) 
            {
                throw new ArgumentException("Weight cannot be negative.");
            }

            double totalWeight = x.Sum(t => t.Item2);
            double cumulativeWeight = 0;
            double cumulativeProb = 0;
            double[] quantiles = new double[probs.Count];
            foreach (Tuple<float, float> nwTuple in x.OrderBy(t => t.Item1))
            {
                cumulativeWeight += nwTuple.Item2;
                cumulativeProb = cumulativeWeight / totalWeight;
                for (int i = 0; i < probs.Count; i++) 
                {
                    if (cumulativeProb <= probs[i]) { quantiles[i] = nwTuple.Item1; }
                }
            }

            return quantiles;
        }

        public static double WeightedMedian(List<Tuple<float, float>> x) 
        {
            return WeightedQuantiles(x, new List<float>() { 0.5f })[0];
        }

        /// <summary>
        /// Weighted interquartile range
        /// </summary>
        public static double WeightedIQR(List<Tuple<float, float>> x)
        {
            double [] iqr = WeightedQuantiles(x, new List<float>() { 0.25f, 0.75f });
            return iqr[1]-iqr[0];
        }


        /// <summary>
        /// Calculates the coefficient of variation (stddev/mean) of a list
        /// Note: Doesn't test if list.Count == 0 or list.Average() == 0
        /// </summary>
        static public double CoefficientOfVariation(List<float> list)
        {
            double average = list.Average();
            double sumOfSquaresOfDifferences = list.Select(val => (val - average) * (val - average)).Sum();
            double sd = Math.Sqrt(sumOfSquaresOfDifferences / list.Count);
            double cv = sd / average;
            return cv;
        }

        /// <summary>
        /// Returns the natural logarithm of an integer's factorial
        /// </summary>
        /// <param name="x">Number to compute the log factorial of.</param>
        /// <returns>The logged factorial.</returns>
        public static double LogFactorial(int x)
        {
            if (CachedLogFactorial == null)
            {
                CachedLogFactorial = new Dictionary<int, double>();
                CachedLogFactorial[0] = 0;
                CachedLogFactorial[1] = 0;
            }
            if (CachedLogFactorial.ContainsKey(x)) return CachedLogFactorial[x];

            // Calculate the log of the gamma function for x + 1. They're equivalent.
            x++;

            double[] c = new double[8] { 0.0833333333, -0.0027777778, 0.0007936508, -0.0005952381, 0.0008417508, -0.0019175269, 0.0064102564, -0.0295506536 };

            const double halfLogTwoPi = 0.91893853320467274178032973640562;

            double z = 1.0 / (x * x);

            double sum = c[7];

            for (int i = 6; i >= 0; i--)
            {
                sum *= z;
                sum += c[i];
            }

            double series = sum / x;

            double logGamma = (x - 0.5) * Math.Log(x) - x + halfLogTwoPi + series;

            CachedLogFactorial[x - 1] = logGamma;

            return logGamma;

        }

        /// <summary>
        /// Returns the 2-norm of vector
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static double TwoNorm(IEnumerable<double> v)
        {
            double norm2 = 0;
            foreach (double x in v)
            {
                norm2 += x * x;
            }

            return Math.Sqrt(norm2);
        }

        /// <summary>
        /// Returns the normalized vector by 2-norm
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static double[] NormalizeBy2Norm(double[] v)
        {
            double[] normalized = new double[v.Length];
            double size = TwoNorm(v);
            if (size == 0)
            {
                Array.Copy(v, normalized, v.Length);
            }
            else
            {
                for (int i = 0; i < v.Length; i++)
                {
                    normalized[i] = v[i] / size;
                }
            }

            return normalized;
        }

        public static double DotProduct(double[] v1, double[] v2)
        {
            if (v1.Length != v2.Length)
            {
                throw new ApplicationException("Vectors must be of the same dimension to calculate dot product.");
            }

            double dotProduct = 0;
            for (int i = 0; i < v1.Length; i++)
            {
                dotProduct += v1[i] * v2[i];
            }

            return dotProduct;
        }

        public static bool AreOrthogonal(double[] v1, double[] v2, double tolerance = 1e-4)
        {
            double dotProduct = DotProduct(v1, v2);
            if (Math.Abs(dotProduct) > tolerance)
                return false;

            return true;
        }

        /// <summary>
        /// Projects v onto axis. Assumes that axis is a unit vector.
        /// </summary>
        /// <param name="v"></param>
        /// <param name="axis"></param>
        /// <returns></returns>
        public static double[] Project(double[] v, double[] axis)
        {
            if (v.Length != axis.Length)
            {
                throw new ApplicationException("Vector and the axis must be of the same dimension.");
            }

            double size = DotProduct(v, axis); // size of the projected vector
            double[] projected = new double[v.Length];
            for (int i = 0; i < v.Length; i++)
            {
                projected[i] = size * axis[i];
            }

            return projected;
        }

        /// <summary>
        /// Projects v onto the space spanned by axes. Assumes that the axes are orthonormal to each other.
        /// </summary>
        /// <param name="v"></param>
        /// <param name="axes"></param>
        /// <returns></returns>
        public static double[] Project(double[] v, IEnumerable<double[]> axes)
        {
            if (axes == null  || !axes.Any())
            {
                throw new ApplicationException("No axes to project onto.");
            }

            if (axes.Any(axis => axis.Length != v.Length))
            {
                throw new ApplicationException("Vector and the axes must be of the same dimension.");
            }

            double[] projected = Project(v, axes.First());
            foreach (double[] axis in axes.Skip(1))
            {
                double[] tmp = Project(v, axis);
                for (int i = 0; i < projected.Length; i++)
                {
                    projected[i] += tmp[i];
                }
            }

            return projected;
        }

        /// <summary>
        /// Projects v to the orthogonal complement of the space spanned by axes. Assumes that the axes are orthonormal to each other.
        /// </summary>
        /// <param name="v"></param>
        /// <param name="axes"></param>
        /// <returns></returns>
        public static double[] ProjectToComplement(double[] v, IEnumerable<double[]> axes)
        {
            double[] projected = Project(v, axes);

            for (int i = 0; i < v.Length; i++)
            {
                projected[i] = v[i] - projected[i];
            }

            return projected;
        }

        public static Dictionary<string, List<GenomicBin>> LoadBedFile(string bedPath, int? gcIndex = null)
        {
            Dictionary<string, List<GenomicBin>> excludedIntervals = new Dictionary<string, List<GenomicBin>>();
            int count = 0;
            using (StreamReader reader = new StreamReader(bedPath))
            {
                while (true)
                {
                    string fileLine = reader.ReadLine();
                    if (fileLine == null) break;
                    string[] bits = fileLine.Split('\t');
                    string chr = bits[0];
                    if (!excludedIntervals.ContainsKey(chr)) excludedIntervals[chr] = new List<GenomicBin>();
                    GenomicBin interval = new GenomicBin();
                    interval.Chromosome = chr;
                    interval.Start = int.Parse(bits[1]);
                    interval.Stop = int.Parse(bits[2]);
                    if (interval.Start < 0)
                    {
                        throw new ApplicationException(String.Format("Start must be non-negative in a BED file: {0}", fileLine));
                    }
                    if (interval.Start >= interval.Stop) // Do not allow empty intervals
                    {
                        throw new ApplicationException(String.Format("Start must be less than Stop in a BED file: {0}", fileLine));
                    }
                    if (gcIndex.HasValue && gcIndex.Value < bits.Length)
                    {
                        interval.GC = int.Parse(bits[gcIndex.Value]);
                    }
                    excludedIntervals[chr].Add(interval);
                    count++;
                }
            }
            Console.WriteLine(">>> Loaded {0} intervals for {1} sequences", count, excludedIntervals.Keys.Count);
            return excludedIntervals;
        }

        static public void SortAndOverlapCheck(Dictionary<string, List<GenomicBin>> intervals, string bedPath)
        {
            List<string> chrs = new List<string>(intervals.Keys);
            for (int chrIndex = 0; chrIndex < chrs.Count; chrIndex++)
            {
                string chr = chrs[chrIndex];
                if (intervals[chr].Count < 2)
                    continue;
                List<GenomicBin> intervalByChrSorted = intervals[chr].OrderBy(o => o.Start).ToList();
                for (int index = 0; index < intervalByChrSorted.Count - 1; index++)
                {
                    if (intervalByChrSorted[index + 1].Start < intervalByChrSorted[index].Stop)
                    {
                        throw new Exception(string.Format($"Bed interval file {bedPath} could only contain non-overlapping intervals"));
                    }
                }
                intervals[chr] = intervalByChrSorted;
            }
        }

        /// <summary>
        /// Estimate the minor allele frequency (MAF) for diploid data at a given coverage level.  This lets us model the fact
        /// that the minor allele frequency will be close to 0.5 but a little bit less (the lower the coverage, the further from 
        /// 0.5 it goes) for ordinary copy-number-2 regions.  
        /// </summary>
        static public double EstimateDiploidMAF(int CopyNumber, float MeanCoverage)
        {
            double expectedCoverageCN1 = MeanCoverage / 2.0f; // Rough guess: Our tumor genome has copy number 2.5
            double expectedCoverage = CopyNumber * expectedCoverageCN1;
            // Empirical model of minor allele frequency by copy number; the two parameters were optimized using 
            // a simple Nelder-Mead procedure, relative to frequencies generated by random sampling
            return 0.5 - 1 / (3.352 * Math.Pow(expectedCoverage, 0.4747));
        }

        static public void PruneVariantFrequencies(List<CanvasSegment> segments, string tempFolder, ref int MinimumVariantFrequenciesForInformativeSegment)
        {
            string debugPath = Path.Combine(tempFolder, "BAFDistributionPerSegment.txt");
            using (StreamWriter writer = new StreamWriter(debugPath))
            {
                foreach (CanvasSegment segment in segments)
                {
                    writer.Write(String.Format("{0}\t{1}\t{2}", segment.Chr, segment.Begin, segment.End));
                    float[] mafs = new float[segment.VariantFrequencies.Count];
                    List<int> zeroIndices = new List<int>();
                    List<int> nonZeroIndices = new List<int>();
                    bool isGreaterThan20 = false;
                    for (int i = 0; i < mafs.Length; i++)
                    {
                        float f = segment.VariantFrequencies[i];
                        mafs[i] = f > 0.5 ? 1 - f : f;
                        if (mafs[i] == 0)
                        {
                            zeroIndices.Add(i);
                        }
                        else
                        {
                            nonZeroIndices.Add(i);
                            if (mafs[i] >= 0.2) { isGreaterThan20 = true; }
                        }
                    }
                    // Segment meets one of the two criteria:
                    //  (1) At least 10% of the alleles have positive MAF
                    //  (2) At least one of the alleles have a MAF >= 0.2
                    if ((float)nonZeroIndices.Count / (float)mafs.Length > 0.1 || isGreaterThan20)
                    {
                        segment.VariantFrequencies = nonZeroIndices.Select(i => segment.VariantFrequencies[i]).ToList();
                        var tmpVFs = segment.VariantFrequencies.Where(v => v > 0.1).ToList(); // heuristic to use only the right mode
                        if (tmpVFs.Count > 0) { segment.VariantFrequencies = tmpVFs; }
                        if (mafs.Length >= MinimumVariantFrequenciesForInformativeSegment) // adjust MinimumVariantFrequenciesForInformativeSegment
                        {
                            MinimumVariantFrequenciesForInformativeSegment = Math.Min(MinimumVariantFrequenciesForInformativeSegment, segment.VariantFrequencies.Count);
                        }
                        writer.Write("\tTrue");
                    }
                    else
                    {
                        writer.Write("\tFalse");
                    }
                    writer.WriteLine("\t" + String.Join(",", mafs));
                }
            }
        }

        /// <summary>
        /// Search for argmin f(x) in [a, b]. Assumes that f is unimodal.
        /// https://en.wikipedia.org/wiki/Golden_section_search
        /// </summary>
        /// <param name="f"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="tol"></param>
        /// <returns></returns>
        public static double GoldenSectionSearch(Func<double, double> f, double a, double b, double tol=1E-5)
        {
            const double goldenRatio = 0.618034; // golden ratio - 1: (Math.Sqrt(5) - 1) / 2

            //  a  c  d  b
            double c = b - goldenRatio * (b - a);
            double d = a + goldenRatio * (b - a);
            double fc = f(c);
            double fd = f(d);
            while (Math.Abs(d - c) > tol)
            {
                if (fc < fd)
                {
                    b = d;
                    d = c;
                    fd = fc;
                    c = b - goldenRatio * (b - a);
                    fc = f(c);
                }
                else
                {
                    a = c;
                    c = d;
                    fc = fd;
                    d = a + goldenRatio * (b - a);
                    fd = f(d);
                }
            }

            return (b + a) / 2;
        }
        /// <summary>
        /// Returns true if a character is an upper or lower case G or C.
        /// </summary>
        /// <param name="c">Character to check.</param>
        /// <returns>True if c is a G or a C (case-insensitive).</returns>
        public static bool IsGC(char c)
        {
            switch (c)
            {
                case 'C': return true;
                case 'G': return true;
                case 'c': return true;
                case 'g': return true;
                default: return false;
            }
        }

        /// <summary>
        /// Is a a subset of b?
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static bool IsSubset<T>(IEnumerable<T> a, IEnumerable<T> b)
        {
            return !a.Except(b).Any();
        }
    }

}
