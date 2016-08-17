using System;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.Linq;
using System.Threading.Tasks;

using Isas.Shared;
using CanvasCommon;

namespace CanvasSmooth
{
    /// <summary>
    /// Smoothing by repeated median filter.
    /// </summary>
    public class CanvasSmooth
    {
        private readonly uint MaxHalfWindowSize;

        public CanvasSmooth(uint maxHalfWindowSize)
        {
            MaxHalfWindowSize = maxHalfWindowSize;
        }

        public int Run(IFileLocation inputFile, IFileLocation outputFile)
        {
            // read input bins
            OrderedDictionary binsByChrom = CanvasIO.GetGenomicBinsByChrom(inputFile.FullName);

            // smooth bins on each chromosome
            RepeatedMedianSmoother smoother = new RepeatedMedianSmoother(MaxHalfWindowSize);
            var chromosomes = binsByChrom.Keys.Cast<string>();
            Dictionary<string, List<GenomicBin>> smoothedBinsByChrom = new Dictionary<string, List<GenomicBin>>();
            Console.WriteLine("Launch smoothing jobs...");
            Console.Out.WriteLine();
            Parallel.ForEach(chromosomes, chrom =>
            {
                smoothedBinsByChrom[chrom] = smoother.Smooth((List<GenomicBin>)binsByChrom[chrom]);
            });
            Console.WriteLine("Completed smoothing jobs.");
            Console.Out.WriteLine();

            // write smoothed bins
            CanvasIO.WriteToTextFile(outputFile.FullName, chromosomes.SelectMany(chrom => smoothedBinsByChrom[chrom]));

            return 0;
        }
    }

    public class RepeatedMedianSmoother
    {
        private readonly uint MaxHalfWindowSize;

        public RepeatedMedianSmoother(uint maxHalfWindowSize)
        {
            MaxHalfWindowSize = maxHalfWindowSize;
        }

        public List<GenomicBin> Smooth(List<GenomicBin> bins)
        {
            IEnumerable<float> counts = bins.Select(b => b.Count);
            IEnumerable<float> smoothedCounts = RepeatedMedianFilter(counts, MaxHalfWindowSize);
            List<GenomicBin> smoothedBins = new List<GenomicBin>();
            smoothedBins.AddRange(Enumerable.Zip(bins, smoothedCounts, (bin, count)
                => new GenomicBin(bin.Chromosome, bin.Start, bin.Stop, bin.GC, count)));
            return smoothedBins;
        }

        private static IEnumerable<float> RepeatedMedianFilter(IEnumerable<float> counts, uint maxHalfWindowSize)
        {
            IEnumerable<float> smoothedCounts = counts;

            for (uint halfWindowSize = 1; halfWindowSize <= maxHalfWindowSize; halfWindowSize++)
            {
                smoothedCounts = CanvasCommon.Utilities.MedianFilter(counts, halfWindowSize);
                counts = smoothedCounts;
            }

            return smoothedCounts;
        }
    }
}
