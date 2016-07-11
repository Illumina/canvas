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

        public IEnumerable<GenomicBin> Run(IFileLocation sampleBedFile, IFileLocation referenceBedFile)
        {
            if (!sampleBedFile.Exists)
                throw new FileNotFoundException(sampleBedFile.FullName + " does not exist.");
            if (!referenceBedFile.Exists)
                throw new FileNotFoundException(referenceBedFile.FullName + " does not exist.");

            double sampleMedian = (new BinCounts(sampleBedFile.FullName, manifest: _manifest)).OnTargetMedianBinCount;
            double referenceMedian = (new BinCounts(referenceBedFile.FullName, manifest: _manifest)).OnTargetMedianBinCount;
            double librarySizeFactor = (sampleMedian > 0 && referenceMedian > 0) ? referenceMedian / sampleMedian : 1;

            using (GzipReader sampleReader = new GzipReader(sampleBedFile.FullName))
            using (GzipReader referenceReader = new GzipReader(referenceBedFile.FullName))
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
                    ratio = sampleCount / referenceCount * librarySizeFactor;
                    yield return new GenomicBin(chrom, start, end, -1, (float)ratio);
                }
            }
        }
    }
}
