using System.Collections.Generic;
using System.Linq;
using CanvasCommon;
using Illumina.Common;
using Isas.SequencingFiles.Bed;

namespace CanvasPedigreeCaller.Visualization
{
    public abstract class BaseNormalizedCoverageCalculator
    {
        /// <summary>
        /// Returns BedGraphEntries for the given CanvasSegments. Can be used with a pre-computed normalization factor (for example,
        /// if the normalization factor was calculated based off of a superset or subset of the segments at hand), or can 
        /// generate the normalization factor on the fly.
        /// </summary>
        /// <param name="segments">List of segments to return normalized coverage bedgraph entries for.</param>
        /// <param name="normalizationFactor">If provided, pre-computed normalization factor to be used for normalizing the segments. If not provided, normalization factor will be calculated based on the provided segments.</param>
        /// <returns></returns>
        public IEnumerable<BedGraphEntry> Calculate(IReadOnlyList<CanvasSegment> segments, double? normalizationFactor = null)
        {
            if (segments.Empty())
                return Enumerable.Empty<BedGraphEntry>();

            if (normalizationFactor == null)
            {
                // If not passed in, calculate normalization factor on the fly using the given segments.
                normalizationFactor = NormalizationCalculator.ComputeNormalizationFactor(segments);
            }

            return segments.SelectMany(segment => GetBedGraphEntries(segment, normalizationFactor.Value));
        }

        protected abstract IEnumerable<BedGraphEntry> GetBedGraphEntries(CanvasSegment segment,
            double normalizationFactor);
    }
}