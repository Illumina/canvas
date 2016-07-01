using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using SequencingFiles;
using Isas.Shared;
using CanvasCommon;

namespace CanvasNormalize
{
    /// <summary>
    /// 
    /// </summary>
    public class LSNormRatioCalculator : IRatioCalculator
    {
        private NexteraManifest _manifest;
        public LSNormRatioCalculator(NexteraManifest manifest)
        {
            _manifest = manifest;
        }

        public void Run(IFileLocation sampleBedFile, IFileLocation referenceBedFile, IFileLocation referencePloidyBedFile,
            IFileLocation outputBedFile)
        {
            if (!sampleBedFile.Exists)
                throw new FileNotFoundException(sampleBedFile.FullName + " does not exist.");
            if (!referenceBedFile.Exists)
                throw new FileNotFoundException(referenceBedFile.FullName + " does not exist.");

            PloidyInfo referencePloidy = null;
            if (referencePloidyBedFile != null && referencePloidyBedFile.Exists)
                referencePloidy = PloidyInfo.LoadPloidyFromBedFile(referencePloidyBedFile.FullName);
            double sampleMedian = (new BinCounts(sampleBedFile.FullName, manifest: _manifest)).OnTargetMedianBinCount;
            double referenceMedian = (new BinCounts(referenceBedFile.FullName, manifest: _manifest)).OnTargetMedianBinCount;
            double librarySizeFactor = (sampleMedian > 0 && referenceMedian > 0) ? referenceMedian / sampleMedian : 1;

            using (GzipReader sampleReader = new GzipReader(sampleBedFile.FullName))
            using (GzipReader referenceReader = new GzipReader(referenceBedFile.FullName))
            using (GzipWriter writer = new GzipWriter(outputBedFile.FullName))
            {
                string referenceLine;
                string sampleLine;
                string[] referenceToks;
                string[] sampleToks;
                double referenceCount;
                double sampleCount;
                double ratio;
                while ((referenceLine = referenceReader.ReadLine()) != null)
                {
                    sampleLine = sampleReader.ReadLine();
                    referenceToks = referenceLine.Split('\t');
                    sampleToks = sampleLine.Split('\t');
                    referenceCount = double.Parse(referenceToks[3]);
                    sampleCount = double.Parse(sampleToks[3]);
                    // The weighted average count of a bin could be less than 1.
                    // Using these small counts for coverage normalization creates large ratios.
                    // It would be better to just drop these bins so we don't introduce too much noise into segmentation and CNV calling.
                    if (referenceCount < 1) { continue; } // skip the bin
                    string chrom = referenceToks[0];
                    int start = int.Parse(referenceToks[1]);
                    int end = int.Parse(referenceToks[2]);
                    // get the normal ploidy from intervalsWithPloidyByChrom
                    double factor = CanvasNormalizeUtilities.CanvasDiploidBinRatioFactor 
                        * CanvasNormalizeUtilities.GetPloidy(referencePloidy, chrom, start, end) / 2.0;
                    ratio = sampleCount / referenceCount * factor * librarySizeFactor;
                    referenceToks[3] = String.Format("{0}", ratio);
                    writer.WriteLine(String.Join("\t", referenceToks));
                }
            }
        }
    }
}
