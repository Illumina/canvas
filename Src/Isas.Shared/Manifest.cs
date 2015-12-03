using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using SequencingFiles;
using ProtoBuf;

namespace Isas.Shared
{
    public static class NexteraManifestUtils
    {
        /// <summary>
        /// Output bed file of targets. These are the intervals after removing the probes
        /// Note that the BED format uses:
        /// 0-based start position (inclusive) and 1-based end position (inclusive)
        /// which is equivalent to saying:
        /// 0-based start position (inclusive) and 0-based end position (exclusive)
        /// </summary>
        public static void WriteTargetBed(NexteraManifest manifest, string outputPath, GenomeMetadata genome)
        {
            using (BgzipOrStreamWriter writer = new BgzipOrStreamWriter(outputPath))
            {
                WriteTargetBed(manifest, writer, genome);
            }
        }

        public static void WriteTargetBed(NexteraManifest manifest, BgzipOrStreamWriter writer, GenomeMetadata genome)
        {
            List<NexteraManifest.ManifestRegion> tempRegions = manifest.Regions;
            if (genome != null)
            {
                tempRegions = new List<NexteraManifest.ManifestRegion>(manifest.Regions);
                Dictionary<string, int> chromsomeIndexLookup = new Dictionary<string, int>();
                //generate chromsome index lookup and sort
                for (int chromosomeIndex = 0; chromosomeIndex < genome.Sequences.Count; chromosomeIndex++)
                {
                    GenomeMetadata.SequenceMetadata sequence = genome.Sequences[chromosomeIndex];
                    chromsomeIndexLookup[sequence.Name] = chromosomeIndex;
                }
                tempRegions.Sort((a, b) => a.CompareTo(b, chromsomeIndexLookup));
            }

                foreach (NexteraManifest.ManifestRegion region in tempRegions)
                {
                TargetInterval interval = region.GetTargetInterval();
                    writer.WriteLine(string.Join("\t", new[]
                    {
                        interval.ReferenceName, 
                        (interval.Begin - 1).ToString(CultureInfo.InvariantCulture), 
                        interval.End.ToString(CultureInfo.InvariantCulture),
                    region.Name //region name is needed for PUMA metrics outputs to generate .coverage.csv file
                    }));
                }
            }

        /// <summary>
        /// Output bed file of regions. Each region spans both probes and the target interval
        /// Note that the BED format uses:
        /// 0-based start position (inclusive) and 1-based end position (inclusive)
        /// which is equivalent to saying:
        /// 0-based start position (inclusive) and 0-based end position (exclusive)
        /// </summary>
        public static void WriteRegionBed(NexteraManifest manifest, string outputPath, GenomeMetadata genome)
        {
            using (BgzipOrStreamWriter writer = new BgzipOrStreamWriter(outputPath))
            {
                WriteRegionBed(manifest, writer, genome);
            }
        }

        public static void WriteRegionBed(NexteraManifest manifest, BgzipOrStreamWriter writer, GenomeMetadata genome)
        {
            List<NexteraManifest.ManifestRegion> tempRegions = manifest.Regions;
            if (genome != null)
            {
                tempRegions = new List<NexteraManifest.ManifestRegion>(manifest.Regions);
                Dictionary<string, int> chromsomeIndexLookup = new Dictionary<string, int>();
                //generate chromsome index lookup and sort
                for (int chromosomeIndex = 0; chromosomeIndex < genome.Sequences.Count; chromosomeIndex++)
                {
                    GenomeMetadata.SequenceMetadata sequence = genome.Sequences[chromosomeIndex];
                    chromsomeIndexLookup[sequence.Name] = chromosomeIndex;
                }
                tempRegions.Sort((a, b) => a.CompareTo(b, chromsomeIndexLookup));
            }

                foreach (NexteraManifest.ManifestRegion region in tempRegions)
                {
                    writer.WriteLine(string.Format("{0}\t{1}\t{2}", region.Chromosome, region.Start - 1, region.End));
                }
            }

        public static void WriteNexteraManifests(NexteraManifest manifest, string path)
        {
            using (StreamWriter writer = new StreamWriter(path))
            {
                WriteNexteraManifests(manifest, writer);
            }
        }

        public static void WriteNexteraManifests(NexteraManifest manifest, TextWriter writer)
        {
                writer.WriteLine("#{0}: {1}", "Manifest Type", "Regions");
            writer.WriteLine("#{0}: {1}", "Target Region Count", manifest.Regions.Count);
                writer.WriteLine("#{0}: {1}", "Date", DateTime.Now.ToShortDateString());
                writer.WriteLine("[Header]");
                //writer.WriteLine("Manifest Version\t1.0");
                if (!string.IsNullOrEmpty(manifest.GenomeName))
                {
                    writer.WriteLine("ReferenceGenome\t" + manifest.GenomeName);
                }
                writer.WriteLine("[Regions]");
                List<string> headers = new List<string>();
                if (manifest.ColumnNames != null && manifest.ColumnNames.Length > 0)
                {
                    foreach (int columnNumber in manifest.ColumnNumbers)
                    {
                        if (columnNumber >= 0 && columnNumber < manifest.ColumnNames.Length)
                        {
                            headers.Add(manifest.ColumnNames[columnNumber]);
                        }
                    }
                    writer.WriteLine(string.Join("\t", headers.ToArray()));
                }

                if (manifest.Regions != null && manifest.Regions.Count > 0)
                {
                    foreach (NexteraManifest.ManifestRegion region in manifest.Regions)
                    {
                        NexteraManifest.ManifestRegion tmpRegion = new NexteraManifest.ManifestRegion(region);

                        List<string> line = new List<string>();
                        foreach (string header in headers)
                        {
                            switch (header.ToLowerInvariant())
                            {
                                case "name":
                                    line.Add(tmpRegion.Name);
                                    break;
                                case "chromosome":
                                    line.Add(tmpRegion.Chromosome);
                                    break;
                                case "start":
                                case "amplicon start":
                                    line.Add(tmpRegion.Start.ToString());
                                    break;
                                case "end":
                                case "amplicon end":
                                    line.Add(tmpRegion.End.ToString());
                                    break;
                                case "startprobelength":
                                case "upstream probe length":
                                    line.Add(tmpRegion.StartProbeLength.ToString());
                                    break;
                                case "endprobelength":
                                case "downstream probe length":
                                    line.Add(tmpRegion.EndProbeLength.ToString());
                                    break;
                                case "groupname":
                                case "group name":
                                case "group":
                                case "ip group":
                                    line.Add(tmpRegion.GroupName);
                                    break;
                            }
                        }
                        writer.WriteLine(string.Join("\t", line.ToArray()));
                    }
                }
            }

        public static NexteraManifest GetUpdatedNexteraManifestsWithNewRegions(NexteraManifest manifest, List<RegionStatistics> regionStats)
        {
            NexteraManifest nexteraManifestsWithNewRegions = new NexteraManifest(manifest);

            // create a dictionary for new regions
            Dictionary<string, RegionStatistics> regionStatsLookup = new Dictionary<string, RegionStatistics>();
            foreach (RegionStatistics regionStat in regionStats)
            {
                regionStatsLookup.Add(regionStat.RegionName, regionStat);
        }

            //update the regions
            if (nexteraManifestsWithNewRegions.Regions != null && nexteraManifestsWithNewRegions.Regions.Count > 0)
            {
                var newRegions = new List<NexteraManifest.ManifestRegion>();
                foreach (NexteraManifest.ManifestRegion region in nexteraManifestsWithNewRegions.Regions)
                {
                    NexteraManifest.ManifestRegion tmpRegion = new NexteraManifest.ManifestRegion(region);
                    if (!regionStatsLookup.ContainsKey(region.Name))
                        continue;

                    tmpRegion.Start = regionStatsLookup[region.Name].StartPosition;
                    tmpRegion.End = regionStatsLookup[region.Name].EndPosition;
                    newRegions.Add(tmpRegion);

                }
                nexteraManifestsWithNewRegions.Regions = newRegions;
            }

            return nexteraManifestsWithNewRegions;
        }
    }

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

    // Note: 2015-10-23: (ImplicitFields = ImplicitFields.AllFields, AsReferenceDefault = true)
    //    causes "ProtoBuf.ProtoException: Possible recursion detected"
    // See TSAW-203
    [ProtoContract(SkipConstructor = true)]
	public class ProbeSetTarget
    {
        #region Members
        [ProtoMember(1)]
        public string AmpliconID;
        [ProtoMember(2)]
        public string AssayID;
        [ProtoMember(3)]
        public int ConcatenatedOffset;
        [ProtoMember(4)]
        public string FullAmpliconSequence;
        [ProtoMember(5)]
        public string FullReverseSequence;
        [ProtoMember(6)]
        public string GeneName;
        [ProtoMember(7, AsReference = true)]
        public AmpliconManifest Manifest;
        [ProtoMember(8)]
        public string Name;
        [ProtoMember(9)]
        public CigarAlignment ReferenceAlignment; // optional: used for targets which recreate an indel relative to the reference
        [ProtoMember(10)]
        public int ReferenceIndex; // used for BAM I/O
        [ProtoMember(11)]
        public string Chromosome { get; set; }
        [ProtoMember(12)]
        public int StartPosition { get; set; } // 1-based
        [ProtoMember(13)]
        public int EndPosition { get; set; } // 1-based, inclusive!
        [ProtoMember(14)]
        public bool ReverseStrand { get; set; }
        [ProtoMember(15)]
        public string Sequence { get; set; }
        [ProtoMember(16)]
        public string ReverseSequence { get; set; }
        [ProtoMember(17, AsReference = true)]
        public ProbeSet ProbeA { get; set; }
        [ProtoMember(18, AsReference = true)]
        public ProbeSet ProbeB { get; set; }
        [ProtoMember(19)]
        public int Index { get; set; }
        [ProtoMember(20)]
        public string Build { get; set; }
        [ProtoMember(21)]
        public string Species { get; set; }
        [ProtoMember(22)]
        public string SoftClipReadDirection { get; set; }
        [ProtoMember(23)]
        public int SoftClipUpstream { get; set; }
        [ProtoMember(24)]
        public int SoftClipDownstream { get; set; }
        #endregion
    }

    // Note: 2015-10-23: (ImplicitFields = ImplicitFields.AllFields, AsReferenceDefault = true)
    //    causes "ProtoBuf.ProtoException: Possible recursion detected"
    // See TSAW-203
    [ProtoContract(SkipConstructor = true)]
    public class ProbeSet
    {
        #region Members
        [ProtoMember(1, AsReference = true)]
        public List<ProbeSetTarget> Targets = new List<ProbeSetTarget>();
        [ProtoMember(2)]
        public string LocusID { get; set; } // assay ID in targeted RNA-Seq
        [ProtoMember(3)]
        public string Species { get; set; }
        [ProtoMember(4)]
        public string BuildID { get; set; }
        [ProtoMember(5)]
        public string SubmittedSequence { get; set; }
        [ProtoMember(6)]
        public string Chromosome { get; set; }
        [ProtoMember(7)]
        public int StartPosition { get; set; }
        [ProtoMember(8)]
        public int EndPosition { get; set; }
        [ProtoMember(9)]
        public bool ReverseStrand { get; set; }
        [ProtoMember(10)]
        public bool SubmittedSequenceReverseStrand { get; set; }
        [ProtoMember(11)]
        public string Read1Tag { get; set; } // upstream 
        [ProtoMember(12)]
        public string Read2Tag { get; set; } // downstream
        #endregion
    }

    /// <summary>
    ///     TargetedRegions keeps track of all of the regions of interest targeted by a given manifest.
    /// </summary>
    public class TargetedRegions
    {
        public Dictionary<string, int[]> CoverageData = new Dictionary<string, int[]>();
        // A, C, G, T, Total, ErrorCount

        public Dictionary<string, GenomicInterval> RegionsByChromosome = new Dictionary<string, GenomicInterval>();
    }

    /// <summary>
    ///     Store the information in the [Header] section
    /// </summary>
    [ProtoContract]
    public class HeaderSection
    {
        [ProtoMember(1)]
        public string Manifest; // Manifest name from DesignStudio
    }
}
