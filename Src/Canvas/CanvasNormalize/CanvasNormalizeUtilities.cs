using System.Collections.Generic;
using System.IO;
using CanvasCommon;
using Illumina.Common.CSV;
using Illumina.Common.FileSystem;
using Isas.SequencingFiles;

namespace CanvasNormalize
{
    public class CanvasNormalizeUtilities
    {
        public static readonly double CanvasDiploidBinRatioFactor = 40;

        public static int GetPloidy(PloidyInfo referencePloidy, string chrom, int start, int end, int defaultPloidy = 2)
        {
            if (referencePloidy == null) { return defaultPloidy; }

            CanvasSegment segment = new CanvasSegment(chrom, start, end, new List<SampleGenomicBin>());

            return referencePloidy.GetReferenceCopyNumber(segment);
        }

        private static IEnumerable<SampleGenomicBin> RatiosToCounts(IEnumerable<SampleGenomicBin> ratios, PloidyInfo referencePloidy)
        {
            foreach (SampleGenomicBin ratio in ratios)
            {
                // get the normal ploidy
                double factor = CanvasDiploidBinRatioFactor * GetPloidy(referencePloidy, ratio.GenomicBin.Chromosome, ratio.Start, ratio.Stop) / 2.0;
                double count = ratio.Count * factor;
                yield return new SampleGenomicBin(ratio.GenomicBin.Chromosome, ratio.Start, ratio.Stop, ratio.GenomicBin.GC, (float)count);
            }
        }

        public static void RatiosToCounts(IEnumerable<SampleGenomicBin> ratios, IFileLocation referencePloidyBedFile,
            IFileLocation outputPath)
        {
            PloidyInfo referencePloidy = null;
            if (referencePloidyBedFile != null && referencePloidyBedFile.Exists)
                referencePloidy = PloidyInfo.LoadPloidyFromBedFile(referencePloidyBedFile.FullName);

            CanvasIO.WriteToTextFile(outputPath.FullName, RatiosToCounts(ratios, referencePloidy));
        }

        /// <summary>
        /// Writes copy-number data (cnd) file.
        /// </summary>
        /// <param name="fragmentCountFile"></param>
        /// <param name="referenceCountFile"></param>
        /// <param name="ratios"></param>
        /// <param name="outputFile"></param>
        public static void WriteCndFile(IFileLocation fragmentCountFile, IFileLocation referenceCountFile,
            IEnumerable<SampleGenomicBin> ratios, IFileLocation outputFile)
        {
            IEnumerable<SampleGenomicBin> fragmentCounts = CanvasIO.IterateThroughTextFile(fragmentCountFile.FullName);
            IEnumerable<SampleGenomicBin> referenceCounts = CanvasIO.IterateThroughTextFile(referenceCountFile.FullName);

            using (var eFragment = fragmentCounts.GetEnumerator())
            using (var eReference = referenceCounts.GetEnumerator())
            using (var eRatio = ratios.GetEnumerator())
            using (FileStream stream = new FileStream(outputFile.FullName, FileMode.Create, FileAccess.Write))
            using (StreamWriter writer = new StreamWriter(stream))
            {
                writer.WriteLine(CSVWriter.GetLine("Fragment Count", "Reference Count", "Chromosome",
                    "Start", "End", "Unsmoothed Log Ratio"));
                while (eFragment.MoveNext() && eReference.MoveNext() && eRatio.MoveNext())
                {
                    // Some bins could have been skipped when calculating the ratios
                    while (!eFragment.Current.IsSameBin(eRatio.Current))
                    {
                        if (!eFragment.MoveNext()) // Ran out of fragment bins
                            throw new Illumina.Common.IlluminaException("Fragment bins and ratio bins are not in the same order.");
                    }
                    while (!eReference.Current.IsSameBin(eRatio.Current))
                    {
                        if (!eReference.MoveNext()) // Ran out of reference bins
                            throw new Illumina.Common.IlluminaException("Reference bins and ratio bins are not in the same order.");
                    }
                    if (!eFragment.Current.IsSameBin(eReference.Current)
                        || !eFragment.Current.IsSameBin(eRatio.Current))
                    {
                        throw new Illumina.Common.IlluminaException("Bins are not in the same order.");
                    }
                    writer.WriteLine(CSVWriter.GetLine(eFragment.Current.Count.ToString(),
                        eReference.Current.Count.ToString(), eFragment.Current.GenomicBin.Chromosome,
                        eFragment.Current.Start.ToString(), eFragment.Current.Stop.ToString(),
                        eRatio.Current.Count.ToString()));
                }
            }
        }
    }
}
