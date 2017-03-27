using System.Collections.Generic;
using System.IO;
using CanvasCommon;
using Illumina.Common.FileSystem;
using Isas.Manifests.NexteraManifest;

namespace CanvasNormalize
{
    public class RawRatioCalculator : IRatioCalculator
    {
        private NexteraManifest _manifest;
        private double _minReferenceCount;
        private double _maxReferenceCount;

        public RawRatioCalculator(NexteraManifest manifest, double minReferenceCount = 1,
            double maxReferecneCount = double.PositiveInfinity)
        {
            _manifest = manifest;
            _minReferenceCount = minReferenceCount;
            _maxReferenceCount = maxReferecneCount;
        }

        public IEnumerable<SampleGenomicBin> Run(IFileLocation sampleBedFile, IFileLocation referenceBedFile)
        {
            if (!sampleBedFile.Exists)
                throw new FileNotFoundException(sampleBedFile.FullName + " does not exist.");
            if (!referenceBedFile.Exists)
                throw new FileNotFoundException(referenceBedFile.FullName + " does not exist.");

            var sampleBins = CanvasIO.IterateThroughTextFile(sampleBedFile.FullName);
            var referenceBins = CanvasIO.IterateThroughTextFile(referenceBedFile.FullName);
            using (var eSampleBins = sampleBins.GetEnumerator())
            using (var eReferenceBins = referenceBins.GetEnumerator())
            {
                while (eSampleBins.MoveNext() && eReferenceBins.MoveNext())
                {
                    var sampleBin = eSampleBins.Current;
                    var referenceBin = eReferenceBins.Current;
                    // Bins with extreme reference counts introduce large variance into the ratios.
                    // It would be better to just drop these bins so we don't introduce too much noise into segmentation and CNV calling.
                    if (referenceBin.Count < _minReferenceCount) { continue; } // skip the bin
                    if (referenceBin.Count > _maxReferenceCount) { continue; } // skip the bin
                    double sampleCount = eSampleBins.Current.Count;
                    double ratio = sampleBin.Count / referenceBin.Count;
                    yield return new SampleGenomicBin(sampleBin.GenomicBin.Chromosome, sampleBin.Start, sampleBin.Stop, sampleBin.GenomicBin.GC, (float)ratio);
                }
            }
        }
    }
}
