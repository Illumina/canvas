using System;
using System.Collections.Generic;
using System.Linq;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Isas.Framework.DataTypes;
using Isas.SequencingFiles;
using Isas.SequencingFiles.Vcf;

namespace CanvasCommon
{
    public class ReferencePloidy
    {
        private readonly Dictionary<string, List<(Interval Interval, int ReferencePloidy)>> _regions;

        private ReferencePloidy(Dictionary<string, List<(Interval Interval, int ReferencePloidy)>> regions)
        {
            _regions = regions;
        }

        /// <summary>
        /// Returns a sequence of adjacent <see cref="ReferencePloidyInterval"/>s that span the provided <paramref name="referenceInterval"/>.<para/>
        /// <see cref="ReferencePloidyInterval"/>s are trimmed to only cover the provided <paramref name="referenceInterval"/><para/>
        /// Adjacent <see cref="ReferencePloidyInterval"/>s have different reference ploidy
        /// </summary>
        public IEnumerable<ReferencePloidyInterval> GetReferencePloidyIntervals(ReferenceInterval referenceInterval)
        {
            string chromosome = referenceInterval.Chromosome;
            if (!_regions.TryGetValue(chromosome, out var regions))
            {
                regions = Enumerable.Empty<(Interval, int)>().ToList();
            }

            var remainingInterval = referenceInterval.Interval;
            foreach (var (interval, ploidy) in regions)
            {
                if (!interval.TryGetOverlap(remainingInterval, out var overlapInterval))
                    continue;

                // handle overhang to the left
                if (remainingInterval.OneBasedStart < interval.OneBasedStart)
                {
                    var diploidInterval = new Interval(remainingInterval.OneBasedStart, interval.OneBasedStart - 1);
                    var diploidReferenceInterval = new ReferenceInterval(chromosome, diploidInterval);
                    yield return new ReferencePloidyInterval(diploidReferenceInterval, 2);
                    remainingInterval = new Interval(interval.OneBasedStart, remainingInterval.OneBasedEnd);
                }

                //handle overlap in the middle
                var overlapReferenceInterval = new ReferenceInterval(chromosome, overlapInterval);
                yield return new ReferencePloidyInterval(overlapReferenceInterval, ploidy);

                if (remainingInterval.OneBasedEnd <= interval.OneBasedEnd)
                {
                    yield break;
                }

                remainingInterval = new Interval(interval.OneBasedEnd + 1, remainingInterval.OneBasedEnd);
            }
            // handle remainingInterval which doesn't overlap any ploidy regions
            var remainingReferenceInterval = new ReferenceInterval(chromosome, remainingInterval);
            yield return new ReferencePloidyInterval(remainingReferenceInterval, 2);
        }

        public static ReferencePloidy Load(VcfReader reader, SampleId sampleId)
        {
            var regions = LoadRegions(reader, sampleId);
            return new ReferencePloidy(regions);
        }

        private static Dictionary<string, List<(Interval Interval, int ReferencePloidy)>> LoadRegions(VcfReader vcfReader, SampleId sampleId)
        {
            var sampleGenotypeColumnIndex = vcfReader.Samples.IndexOf(sampleId.ToString());
            if (sampleGenotypeColumnIndex == -1)
                throw new ArgumentException($"VCF does not contain genotype column for sample '{sampleId}'");

            return vcfReader
                .GetVariants()
                .GroupByAdjacent(entry => entry.ReferenceName)
                .Select(kvp => (kvp.Key, LoadReferencePloidyIntervals(kvp.Value, sampleGenotypeColumnIndex)))
                .ToDictionary();
        }

        private static List<(Interval Interval, int ReferencePloidy)> LoadReferencePloidyIntervals(IEnumerable<VcfVariant> entries, int sampleGenotypeColumnIndex)
        {
            var regions = entries.Select(entry => LoadReferencePloidyInterval(entry, sampleGenotypeColumnIndex)).ToList();
            return GetContinuousRegions(regions).Where(region => region.ReferencePloidy != 2).ToList();
        }

        private static IEnumerable<(Interval Interval, int ReferencePloidy)> GetContinuousRegions(List<(Interval Interval, int ReferencePloidy)> regions)
        {
            if (!regions.Any()) yield break;
            var (currentInterval, currentPloidy) = regions.First();
            foreach (var (nextInterval, nextPloidy) in regions.Skip(1))
            {
                if (currentInterval.Overlaps(nextInterval))
                    throw new ArgumentException(
                        $"Error in Ploidy VCF. Found overlapping intervals '{currentInterval}' and '{nextInterval}'");
                if (currentInterval.CompareTo(nextInterval) > 0)
                    throw new ArgumentException(
                        $"Error in Ploidy VCF. Missorted intervals '{currentInterval}' followed by '{nextInterval}'");

                if (currentInterval.LeftAdjacent(nextInterval) && currentPloidy == nextPloidy)
                {
                    currentInterval = new Interval(currentInterval.OneBasedStart, nextInterval.OneBasedEnd);
                    continue;
                }
                yield return (currentInterval, currentPloidy);
                currentInterval = nextInterval;
                currentPloidy = nextPloidy;
            }
            yield return (currentInterval, currentPloidy);
        }

        private static (Interval Interval, int ReferencePloidy) LoadReferencePloidyInterval(VcfVariant entry, int sampleGenotypeColumnIndex)
        {
            var genotypeColumn = entry.GenotypeColumns[sampleGenotypeColumnIndex];
            if (!genotypeColumn.TryGetValue("CN", out var cnValue))
                throw new ArgumentException($"Missing CN field in genotype column. Vcf entry: {entry}");
            if (!uint.TryParse(cnValue, out var referencePloidy))
                throw new ArgumentException($"CN field must be an unsigned integer. Vcf entry: {entry}");

            if (!entry.InfoFields.TryGetValue("END", out var endValue))
                throw new ArgumentException($"Missing END field in INFO column. Vcf entry: {entry}");
            if (!uint.TryParse(endValue, out var end))
                throw new ArgumentException($"END field must be an unsigned integer. Vcf entry: {entry}");

            var start = entry.ReferencePosition;
            if (entry.VariantAlleles[0].StartsWith("<"))
                start += 1; // vcf spec says if ALT is symbolic then POS is before the record (i.e. padding base)
            return (new Interval(start, (int)end), (int)referencePloidy);
        }

        internal static ReferencePloidy Load(IFileLocation vcf, SampleId sampleId)
        {
            using (var reader = new VcfReader(vcf.FullName))
            {
                try
                {
                    return Load(reader, sampleId);
                }
                catch (Exception e)
                {
                    throw new ArgumentException($"Failed to load ploidy VCF at '{vcf}'. Error: {e.Message}", e);
                }
            }
        }
    }
}