using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using CanvasCommon;
using Illumina.Common.FileSystem;

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
            var binsByChrom = CanvasIO.GetGenomicBinsByChrom(inputFile.FullName);

            // smooth bins on each chromosome
            RepeatedMedianSmoother smoother = new RepeatedMedianSmoother(MaxHalfWindowSize);
            var chromosomes = binsByChrom.Keys;
            ConcurrentDictionary<string, List<SampleGenomicBin>> smoothedBinsByChrom = new ConcurrentDictionary<string, List<SampleGenomicBin>>();
            Console.WriteLine("Launch smoothing jobs...");
            Parallel.ForEach(chromosomes, chrom =>
            {
                smoothedBinsByChrom[chrom] = smoother.Smooth(binsByChrom[chrom]);
            });
            Console.WriteLine("Completed smoothing jobs.");

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

        public List<SampleGenomicBin> Smooth(List<SampleGenomicBin> bins)
        {
            IEnumerable<float> counts = bins.Select(b => b.Count);
            IEnumerable<float> smoothedCounts = RepeatedMedianFilter(counts, MaxHalfWindowSize);
            List<SampleGenomicBin> smoothedBins = new List<SampleGenomicBin>();
            smoothedBins.AddRange(Enumerable.Zip(bins, smoothedCounts, (bin, count)
                => new SampleGenomicBin(bin.GenomicBin.Chromosome, bin.Start, bin.Stop, bin.GenomicBin.GC, count)));
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
