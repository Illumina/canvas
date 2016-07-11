using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Isas.Shared;
using SequencingFiles;
using CanvasCommon;

namespace CanvasNormalize
{
    public class CanvasNormalizeUtilities
    {
        public static readonly double CanvasDiploidBinRatioFactor = 40;

        public static int GetPloidy(PloidyInfo referencePloidy, string chrom, int start, int end, int defaultPloidy = 2)
        {
            if (referencePloidy == null) { return defaultPloidy; }

            CanvasSegment segment = new CanvasSegment(chrom, start, end, new List<float>());

            return referencePloidy.GetReferenceCopyNumber(segment);
        }

        public static void RatiosToCounts(IEnumerable<GenomicBin> ratios, IFileLocation referencePloidyBedFile,
            IFileLocation outputPath)
        {
            PloidyInfo referencePloidy = null;
            if (referencePloidyBedFile != null && referencePloidyBedFile.Exists)
                referencePloidy = PloidyInfo.LoadPloidyFromBedFile(referencePloidyBedFile.FullName);

            using (GzipWriter writer = new GzipWriter(outputPath.FullName))
            {
                foreach (GenomicBin ratio in ratios)
                {
                    // get the normal ploidy
                    double factor = CanvasDiploidBinRatioFactor * GetPloidy(referencePloidy, ratio.Chromosome, ratio.Start, ratio.Stop) / 2.0;
                    double count = ratio.Count * factor;
                    writer.WriteLine(String.Join("\t", ratio.Chromosome, ratio.Start, ratio.Stop, count));
                }
            }
        }
    }
}
