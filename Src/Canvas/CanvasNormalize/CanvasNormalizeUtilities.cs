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

        public static void RatiosToCounts(IEnumerable<Tuple<string, int, int, double>> ratios, IFileLocation referencePloidyBedFile,
            IFileLocation outputPath)
        {
            PloidyInfo referencePloidy = null;
            if (referencePloidyBedFile != null && referencePloidyBedFile.Exists)
                referencePloidy = PloidyInfo.LoadPloidyFromBedFile(referencePloidyBedFile.FullName);

            using (GzipWriter writer = new GzipWriter(outputPath.FullName))
            {
                //writer.WriteLine(String.Join("\t", referenceToks));
                foreach (var ratio in ratios)
                {
                    string chrom = ratio.Item1;
                    int start = ratio.Item2;
                    int end = ratio.Item3;
                    // get the normal ploidy
                    double factor = CanvasDiploidBinRatioFactor * GetPloidy(referencePloidy, chrom, start, end) / 2.0;
                    double count = ratio.Item4 * factor;
                    writer.WriteLine(String.Join("\t", chrom, start, end, count));
                }
            }
        }
    }
}
