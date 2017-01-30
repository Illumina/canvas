using System.Collections.Generic;
using System.Linq;
using CanvasCommon;
using Isas.Manifests.NexteraManifest;
using Isas.SequencingFiles;

namespace CanvasNormalize
{
    class BinCounts
    {
        private string BinnedPath;
        private NexteraManifest Manifest;
        private List<double> Counts = null;
        private List<int> OnTargetIndices = null;
        private double? onTargetMedianBinCount = null;

        public List<double> AllCounts
        {
            get { return Counts; }
        }

        public double MedianBinCount
        {
            get { return Utilities.Median(Counts); }
        }

        /// <summary>
        /// Return the on-target bin counts if a manifest is provided.
        /// Otherwise, return all the bin counts.
        /// </summary>
        public List<double> OnTargetCounts
        {
            get
            {
                if (OnTargetIndices == null)
                {
                    return Counts;
                }
                else
                {
                    List<double> onTargetCounts = OnTargetIndices.Select(idx => Counts[idx]).ToList();
                    onTargetMedianBinCount = Utilities.Median(onTargetCounts);
                    return onTargetCounts;
                }
            }
        }

        public double OnTargetMedianBinCount
        {
            get
            {
                if (onTargetMedianBinCount.HasValue)
                {
                    return onTargetMedianBinCount.Value;
                }
                else
                {
                    return Utilities.Median(OnTargetCounts);
                }
            }
        }

        public BinCounts(string binnedPath, NexteraManifest manifest = null)
        {
            BinnedPath = binnedPath;
            Manifest = manifest;
            LoadBinCounts(binnedPath, manifest);
        }

        public BinCounts(IEnumerable<SampleGenomicBin> bins, NexteraManifest manifest = null)
        {
            Manifest = manifest;
            LoadBinCounts(bins, manifest);
        }

        private void LoadBinCounts(string binnedPath, NexteraManifest manifest)
        {
            if (manifest == null)
            {
                LoadBinCounts(binnedPath, out Counts);
            }
            else
            {
                LoadBinCounts(binnedPath, manifest, out Counts, out OnTargetIndices);
            }
        }

        private void LoadBinCounts(IEnumerable<SampleGenomicBin> bins, NexteraManifest manifest)
        {
            if (manifest == null)
            {
                Counts = bins.Select(bin => (double)bin.Count).ToList();
            }
            else
            {
                LoadBinCounts(bins, manifest, out Counts, out OnTargetIndices);
            }
        }

        private static void LoadBinCounts(string binnedPath, out List<double> binCounts)
        {
            binCounts = new List<double>();

            using (GzipReader reader = new GzipReader(binnedPath))
            {
                string line;
                string[] toks;
                while ((line = reader.ReadLine()) != null)
                {
                    toks = line.Split('\t');
                    binCounts.Add(double.Parse(toks[3]));
                }
            }
        }

        private static void LoadBinCounts(IEnumerable<SampleGenomicBin> bins, NexteraManifest manifest,
            out List<double> binCounts, out List<int> onTargetIndices)
        {
            binCounts = new List<double>();
            onTargetIndices = new List<int>();

            var regionsByChrom = manifest.GetManifestRegionsByChromosome();
            string currChrom = null;
            List<NexteraManifest.ManifestRegion> regions = null; // 1-based regions
            int regionIndex = -1;
            bool onTarget = false;
            int binIdx = 0;
            foreach (var bin in bins)
            {
                if (currChrom != bin.GenomicBin.Chromosome)
                {
                    currChrom = bin.GenomicBin.Chromosome;
                    onTarget = false;
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
                    onTarget = true;
                }
                else
                {
                    onTarget = false;
                }

                if (onTarget) { onTargetIndices.Add(binIdx); }

                binCounts.Add(bin.Count);
                binIdx++;
            }
        }

        private static void LoadBinCounts(string binnedPath, NexteraManifest manifest, out List<double> binCounts,
            out List<int> onTargetIndices)
        {
            LoadBinCounts(CanvasIO.IterateThroughTextFile(binnedPath), manifest, out binCounts, out onTargetIndices);
        }
    }
}
