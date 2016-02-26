using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using SequencingFiles;
using CanvasCommon;

namespace CanvasClean
{
    public static class EnrichmentUtilities
    {
        /// <summary>
        /// Get the on-target bins by intersecting the manifest.
        /// </summary>
        /// <param name="bins"></param>
        /// <param name="manifest"></param>
        /// <returns></returns>
        public static IEnumerable<GenomicBin> GetOnTargetBins(IEnumerable<GenomicBin> bins, NexteraManifest manifest)
        {
            var regionsByChrom = manifest.GetManifestRegionsByChromosome();
            string currChrom = null;
            List<NexteraManifest.ManifestRegion> regions = null; // 1-based regions
            int regionIndex = -1;
            bool offTarget = true;
            foreach (GenomicBin bin in bins) // 0-based bins
            {
                if (currChrom != bin.Chromosome)
                {
                    currChrom = bin.Chromosome;
                    offTarget = true;
                    if (!regionsByChrom.ContainsKey(currChrom))
                    {
                        regions = null;
                    }
                    else
                    {
                        regions = regionsByChrom[currChrom];
                        regionIndex = 0;
                    }
                }
                while (regions != null && regionIndex < regions.Count && regions[regionIndex].End < bin.Start + 1)
                {
                    regionIndex++;
                }
                if (regions != null && regionIndex < regions.Count && regions[regionIndex].Start <= bin.Stop) // overlap
                {
                    offTarget = false;
                }
                else
                {
                    offTarget = true;
                }

                if (offTarget) { continue; } // ignore off-target bins

                yield return bin;
            }
        }

        public static readonly int numberOfGCbins = 101;
        /// <summary>
        /// Assumes the bins are sorted by genomic coordinates
        /// </summary>
        /// <param name="bins">Bins whose counts are to be normalized</param>
        /// <param name="countsByGC">An array of lists. Each array element (0-100) will hold a list of counts whose bins have the same GC content.</param>
        /// <param name="counts">Will hold all of the autosomal counts present in 'bins'</param>
        public static void GetCountsByGC(List<GenomicBin> bins, NexteraManifest manifest, out List<float>[] countsByGC, out List<float> counts)
        {
            countsByGC = new List<float>[numberOfGCbins];
            counts = new List<float>(bins.Count);

            // Initialize the lists
            for (int i = 0; i < countsByGC.Length; i++)
                countsByGC[i] = new List<float>();

            foreach (GenomicBin bin in manifest == null ? bins : EnrichmentUtilities.GetOnTargetBins(bins, manifest))
            {
                if (!GenomeMetadata.SequenceMetadata.IsAutosome(bin.Chromosome)) { continue; }

                // Put the observed count in the GC-appropriate list.
                countsByGC[bin.GC].Add(bin.Count);

                // Add to the global list of counts.
                counts.Add(bin.Count);
            }
        }
    }
}
