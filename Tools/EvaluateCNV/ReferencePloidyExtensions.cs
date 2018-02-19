using System;
using System.Linq;
using CanvasCommon;
using Isas.SequencingFiles;

namespace EvaluateCNV
{
    public static class ReferencePloidyExtensions
    {
        /// <summary>
        /// Returns the reference ploidy for <paramref name="referenceInterval"/><para/>
        /// If <paramref name="referenceInterval"/> spans regions with different ploidy an exception is thrown
        /// </summary>
        public static int GetSingleReferencePloidy(this ReferencePloidy referencePloidy, ReferenceInterval referenceInterval)
        {
            var referencePloidies = referencePloidy.GetReferencePloidyIntervals(referenceInterval).ToList();
            if (referencePloidies.Count != 1)
            {
                var differentPloidyIntervals = referencePloidies.Select(ploidy => ploidy.Interval);
                throw new ArgumentException(
                    $"Reference interval '{referenceInterval}' overlaps regions with different ploidy: '{string.Join(", ", differentPloidyIntervals)}'");
            }

            return referencePloidies.Single().ReferencePloidy;
        }
    }
}