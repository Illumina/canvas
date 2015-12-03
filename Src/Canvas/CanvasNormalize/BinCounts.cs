using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using CanvasCommon;
using SequencingFiles;

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
            LoadBinCounts();
        }

        private void LoadBinCounts()
        {
            if (Manifest == null)
            {
                LoadBinCounts(BinnedPath, out Counts);
            }
            else
            {
                LoadBinCounts(BinnedPath, Manifest, out Counts, out OnTargetIndices);
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

        private static void LoadBinCounts(string binnedPath, NexteraManifest manifest, out List<double> binCounts,
            out List<int> onTargetIndices)
        {
            binCounts = new List<double>();
            onTargetIndices = new List<int>();

            var regionsByChrom = manifest.GetManifestRegionsByChromosome();
            string currChrom = null;
            List<NexteraManifest.ManifestRegion> regions = null; // 1-based regions
            int regionIndex = -1;
            bool onTarget = false;
            using (GzipReader reader = new GzipReader(binnedPath))
            {
                string line;
                string[] toks;
                int binIdx = 0;
                while ((line = reader.ReadLine()) != null)
                {
                    toks = line.Split('\t');
                    string chrom = toks[0];
                    int start = int.Parse(toks[1]); // 0-based, inclusive
                    int stop = int.Parse(toks[2]); // 0-based, exclusive
                    if (currChrom != chrom)
                    {
                        currChrom = chrom;
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
                    while (regions != null && regionIndex < regions.Count && regions[regionIndex].End < start + 1)
                    {
                        regionIndex++;
                    }
                    if (regions != null && regionIndex < regions.Count && regions[regionIndex].Start <= stop) // overlap
                    {
                        onTarget = true;
                    }
                    else
                    {
                        onTarget = false;
                    }

                    if (onTarget) { onTargetIndices.Add(binIdx); }

                    binCounts.Add(double.Parse(toks[3]));
                    binIdx++;
                }
            }
        }
    }
}
