using System.Collections.Generic;

namespace Isas.Shared
{
    /// <summary>
    ///     TargetedRegions keeps track of all of the regions of interest targeted by a given manifest.
    /// </summary>
    public class TargetedRegions
    {
        public Dictionary<string, int[]> CoverageData = new Dictionary<string, int[]>();
        // A, C, G, T, Total, ErrorCount

        public Dictionary<string, GenomicInterval> RegionsByChromosome = new Dictionary<string, GenomicInterval>();
    }
}