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

            var sampleBins = CanvasIO.IterateThroughTextFile(sampleBedFile.FullName);
            var referenceBins = CanvasIO.IterateThroughTextFile(referenceBedFile.FullName);
            double sampleMedian = (new BinCounts(sampleBins, manifest: _manifest)).OnTargetMedianBinCount;
            double referenceMedian = (new BinCounts(referenceBins, manifest: _manifest)).OnTargetMedianBinCount;
            double librarySizeFactor = (sampleMedian > 0 && referenceMedian > 0) ? referenceMedian / sampleMedian : 1;

            using (var eSampleBins = sampleBins.GetEnumerator())
            using (var eReferenceBins = referenceBins.GetEnumerator())
            {
                while (eSampleBins.MoveNext() && eReferenceBins.MoveNext())
                {
                    var sampleBin = eSampleBins.Current;
                    var referenceBin = eReferenceBins.Current;
                    // The weighted average count of a bin could be less than 1.
                    // Using these small counts for coverage normalization creates large ratios.
                    // It would be better to just drop these bins so we don't introduce too much noise into segmentation and CNV calling.
                    if (referenceBin.Count < 1) { continue; } // skip the bin
                    double ratio = sampleBin.Count / referenceBin.Count * librarySizeFactor;
                    yield return new GenomicBin(sampleBin.Chromosome, sampleBin.Start, sampleBin.Stop, sampleBin.GC, (float)ratio);
                }
            }
        }
    }
}
