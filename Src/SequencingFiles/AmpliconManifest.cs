using System;
using System.Collections.Generic;
using System.IO;
using ProtoBuf;
using SequencingFiles;

namespace Isas.Shared
{
    /// <summary>
    ///     An amplicon manifest captures information about our amplification probes, the genomic region they target,
    ///     and any other off-target regions we may end up sequencing.  Used for TSCA (the Custom Amplicon workflow)
    /// </summary>
    [Serializable]
    [ProtoContract(SkipConstructor = true, ImplicitFields = ImplicitFields.AllFields, AsReferenceDefault = true)]
    public class AmpliconManifest
    {
        #region Members
        public HeaderSection HeaderSection;
        public string Name; // File name including any extension (normally .txt)
        public Dictionary<string, int> PositionMapping; // Chromosome:Position -> Concatenated position
        public ProbeSet[] Probes;
        public int PseudogenomeLength;
        public ProbeSetTarget[] Targets;
        public Dictionary<string, List<GenomicInterval>> Intervals = new Dictionary<string, List<GenomicInterval>>();
        public string PrettyName { get; private set; }
        #endregion

        // constructor
        public AmpliconManifest(string prettyName)
        {
            PrettyName = prettyName;
            HeaderSection = new HeaderSection();
        }

        public long CalculateTotalRegionLength(bool excludeExpectedOffTarget = true)
        {
            long manifestLength = 0;
            Dictionary<string, int> chromosomeDict = new Dictionary<string, int>();
            List<List<Tuple<long, long>>> targets = new List<List<Tuple<long, long>>>();

            foreach (ProbeSetTarget target in Targets)
            {
                if (excludeExpectedOffTarget && target.Index > 1) continue;
                if (!chromosomeDict.ContainsKey(target.Chromosome))
                {
                    chromosomeDict.Add(target.Chromosome, chromosomeDict.Count);
                    targets.Add(new List<Tuple<long, long>>());
                }
                targets[chromosomeDict[target.Chromosome]].Add(new Tuple<long, long>(target.StartPosition, target.EndPosition));
            }

            foreach (List<Tuple<long, long>> targetsPerChromosome in targets)
            {
                targetsPerChromosome.Sort((x, y) => x.Item1.CompareTo(y.Item1));
                long endPosition = 0;
                foreach (Tuple<long, long> targetStartStop in targetsPerChromosome)
                {
                    long startPosition = targetStartStop.Item1;
                    // make sure to account for overlapping regions
                    if (endPosition >= startPosition)
                        startPosition = endPosition + 1;
                    if (targetStartStop.Item2 > endPosition)
                        endPosition = targetStartStop.Item2;
                    if (endPosition - startPosition > 0)
                        manifestLength += endPosition - startPosition + 1;
                }
            }

            return manifestLength;
        }
    }
}