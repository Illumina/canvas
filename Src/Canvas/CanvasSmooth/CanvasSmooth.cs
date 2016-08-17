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
            List<string> chromosomes = new List<string>();
            Dictionary<string, List<GenomicBin>> smoothedBinsByChrom = new Dictionary<string, List<GenomicBin>>();
            List<SmoothTask> tasks = new List<SmoothTask>();
            foreach (string chrom in binsByChrom.Keys)
            {
                chromosomes.Add(chrom);
                smoothedBinsByChrom[chrom] = new List<GenomicBin>();
                SmoothTask task = new SmoothTask(MaxHalfWindowSize, (List<GenomicBin>)binsByChrom[chrom], smoothedBinsByChrom[chrom]);
                tasks.Add(task);
            }
            Console.WriteLine("Launch smoothing jobs...");
            Console.Out.WriteLine();
            Parallel.ForEach(tasks, t => t.Run());
            Console.WriteLine("Completed smoothing jobs.");
            Console.Out.WriteLine();

            // write smoothed bins
            CanvasIO.WriteToTextFile(outputFile.FullName, chromosomes.SelectMany(chrom => smoothedBinsByChrom[chrom]));

            return 0;
        }

        public class SmoothTask
        {
            private readonly uint MaxHalfWindowSize;
            private readonly List<GenomicBin> Bins;
            private List<GenomicBin> SmoothedBins;

            public SmoothTask(uint maxHalfWindowSize, List<GenomicBin> bins, List<GenomicBin> smoothedBins)
            {
                MaxHalfWindowSize = maxHalfWindowSize;
                Bins = bins;
                SmoothedBins = smoothedBins;
            }

            public void Run()
            {
                IEnumerable<float> counts = Bins.Select(b => b.Count);
                IEnumerable<float> smoothedCounts = RepeatedMedianFilter(counts, MaxHalfWindowSize);
                SmoothedBins.AddRange(Enumerable.Zip(Bins, smoothedCounts, (bin, count) 
                    => new GenomicBin(bin.Chromosome, bin.Start, bin.Stop, bin.GC, count)));
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
}
