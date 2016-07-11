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

        public IEnumerable<GenomicBin> Run(IFileLocation sampleBedFile, IFileLocation referenceBedFile)
        {
            if (!sampleBedFile.Exists)
                throw new FileNotFoundException(sampleBedFile.FullName + " does not exist.");
            if (!referenceBedFile.Exists)
                throw new FileNotFoundException(referenceBedFile.FullName + " does not exist.");

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
                    // Bins with extreme reference counts introduce large variance into the ratios.
                    // It would be better to just drop these bins so we don't introduce too much noise into segmentation and CNV calling.
                    if (referenceCount < _minReferenceCount) { continue; } // skip the bin
                    if (referenceCount > _maxReferenceCount) { continue; } // skip the bin
                    string chrom = referenceToks[0];
                    int start = int.Parse(referenceToks[1]);
                    int end = int.Parse(referenceToks[2]);
                    ratio = sampleCount / referenceCount;
                    yield return new GenomicBin(chrom, start, end, -1, (float)ratio);
                }
            }
        }
    }
}
