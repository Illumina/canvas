using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using SequencingFiles;
using ProtoBuf;

namespace SequencingFiles
{

    /// <summary>
    /// This manifest format is used for the PCR Amplicon and Enrichment workflows (Wave, Neo).  Note that it's
    /// not necessarily used for nextera sample prep data, despite the word "Nextera" in the name :)
    /// </summary>
    public class NexteraManifest
    {
        #region Members

        public readonly string[] ColumnNames = new[]
        {
            "Name", "Chromosome", "Amplicon Start", "Amplicon End",
            "Upstream Probe Length", "Downstream Probe Length", "Group"
        };

        public readonly int[] ColumnNumbers = new int[7];

        private readonly bool[] ColumnRequiredFlag = new[] { true, true, true, true, false, false, false };
        public string GenomeName;
        public List<ManifestRegion> Regions = new List<ManifestRegion>();
        public string Name { get; set; }
        public int? CanvasBinSize = null;
        public string CanvasControlBinnedPath = null;
        public string CanvasSexChromosomeKaryotype = null;
        public bool CanvasControlAvailable { get { return CanvasBinSize.HasValue && (CanvasControlBinnedPath != null) && (CanvasSexChromosomeKaryotype != null); } }
        #endregion

        public delegate void ErrorHandler(string message);

        public NexteraManifest(NexteraManifest existingManifest)
        {
            ColumnNames = existingManifest.ColumnNames.ToArray();
            ColumnNumbers = existingManifest.ColumnNumbers.ToArray();
            ColumnRequiredFlag = existingManifest.ColumnRequiredFlag.ToArray();
            GenomeName = existingManifest.GenomeName;
            Regions = existingManifest.Regions.ToList();
            Name = existingManifest.Name;
            CanvasBinSize = existingManifest.CanvasBinSize;
            CanvasControlBinnedPath = existingManifest.CanvasControlBinnedPath;
        }

        public NexteraManifest(string filePath, HashSet<string> ExcludeRegionGroups, ErrorHandler Error)
        {
            string manifestName = Path.GetFileName(filePath);
            using (Stream inputStream = new FileStream(filePath, FileMode.Open, FileAccess.Read))
            {
                init(inputStream, manifestName, ExcludeRegionGroups, Error);
            }
        }

        public NexteraManifest(Stream inputStream, string manifestName, HashSet<string> ExcludeRegionGroups, ErrorHandler Error)
        {
            init(inputStream, manifestName, ExcludeRegionGroups, Error);
        }

        protected void init(Stream inputStream, string manifestName, HashSet<string> ExcludeRegionGroups, ErrorHandler Error)
        {
            HashSet<string> GroupNamesSeen = new HashSet<string>();
            Name = manifestName;
            int lineNumber = -1;
            bool inHeader = false;
            HashSet<string> regionNames = new HashSet<string>();
            for (int fieldIndex = 0; fieldIndex < ColumnNames.Length; fieldIndex++)
            {
                ColumnNumbers[fieldIndex] = -1; // Marker: -1 means no column header seen
            }
            int maxRequiredColumn = -1;
            using (StreamReader reader = new StreamReader(inputStream))
            {
                while (true)
                {
                    string fileLine = reader.ReadLine();
                    if (fileLine == null) break;
                    fileLine = fileLine.TrimEnd();
                    if (fileLine.Length == 0) continue;
                    if (inHeader)
                    {
                        string[] headerBits = fileLine.Split('\t');
                        if (headerBits[0].ToLowerInvariant() == "referencegenome")
                            GenomeName = headerBits[1];
                    }

                    if (fileLine.ToLowerInvariant() == "[header]")
                    {
                        inHeader = true;
                        continue;
                    }

                    if (fileLine.ToLowerInvariant() == "[regions]")
                    {
                        inHeader = false;
                        lineNumber = 0;
                        continue;
                    }
                    string[] bits = fileLine.Split('\t');
                    if (lineNumber == 0)
                    {
                        // This is the header line!  Figure out which column each field goes to:
                        for (int bitIndex = 0; bitIndex < bits.Length; bitIndex++)
                        {
                            int fieldIndex = -1;
                            switch (bits[bitIndex].ToLowerInvariant())
                            {
                                case "name":
                                    fieldIndex = 0;
                                    break;
                                case "chromosome":
                                    fieldIndex = 1;
                                    break;
                                case "amplicon start":
                                    fieldIndex = 2;
                                    break;
                                case "start":
                                    fieldIndex = 2;
                                    //renamed column name if the start case is present in stead of amplicon start. Not necessary.
                                    ColumnNames[fieldIndex] = "Start";
                                    break;
                                case "amplicon end":
                                    fieldIndex = 3;
                                    break;
                                case "end":
                                    fieldIndex = 3;
                                    //renamed column name if the end case is present instead of amplicon end. Not necessary.
                                    ColumnNames[fieldIndex] = "End";
                                    break;
                                case "upstream probe length":
                                    fieldIndex = 4;
                                    break;
                                case "downstream probe length":
                                    fieldIndex = 5;
                                    break;
                                case "ip group": // just in case!
                                case "group":
                                    fieldIndex = 6;
                                    break;
                                default:
                                    break; // Column names we don't recognize are ignored silently
                            }
                            if (fieldIndex >= 0)
                            {
                                if (ColumnNumbers[fieldIndex] != -1)
                                {
                                    throw new Exception(string.Format("Duplicate column header {0} in manifest {1}",
                                        ColumnNames[bitIndex], manifestName));
                                }
                                ColumnNumbers[fieldIndex] = bitIndex;
                                if (ColumnRequiredFlag[fieldIndex])
                                    maxRequiredColumn = Math.Max(maxRequiredColumn, bitIndex);
                            }
                        }
                        for (int fieldIndex = 0; fieldIndex < ColumnNames.Length; fieldIndex++)
                        {
                            if (ColumnRequiredFlag[fieldIndex] && ColumnNumbers[fieldIndex] < 0)
                            {
                                throw new Exception(string.Format("Missing column header {0} in manifest {1}", ColumnNames[fieldIndex], manifestName));
                            }
                        }
                    }
                    if (lineNumber != -1)
                    {
                        lineNumber++;
                        if (lineNumber > 1)
                        {
                            ManifestRegion region = new ManifestRegion { Name = bits[ColumnNumbers[0]] };
                            if (regionNames.Contains(region.Name))
                            {
                                throw new Exception(string.Format(
                                    "Error: Duplicate region name {0} in manifest {1}.  Region names must be unique",
                                    region.Name, manifestName));
                            }
                            regionNames.Add(region.Name);
                            if (bits.Length <= maxRequiredColumn)
                            {
                                throw new Exception(string.Format(
                                    "Error: Manifest file {0} has an entry '{1}' with incorrect number of columns.  Please correct the manifest file syntax.",
                                    manifestName, fileLine));
                            }
                            region.Chromosome = bits[ColumnNumbers[1]];
                            region.Start = int.Parse(bits[ColumnNumbers[2]]);
                            region.End = int.Parse(bits[ColumnNumbers[3]]);
                            if (ColumnNumbers[4] >= 0 && bits.Length > ColumnNumbers[4])
                            {
                                region.StartProbeLength = int.Parse(bits[ColumnNumbers[4]]);
                            }
                            if (ColumnNumbers[5] >= 0 && bits.Length > ColumnNumbers[5])
                            {
                                region.EndProbeLength = int.Parse(bits[ColumnNumbers[5]]);
                            }
                            if (ColumnNumbers[6] >= 0 && bits.Length > ColumnNumbers[6])
                            {
                                region.GroupName = bits[ColumnNumbers[6]].Trim();
                                GroupNamesSeen.Add(region.GroupName);
                                if (ExcludeRegionGroups != null && ExcludeRegionGroups.Contains(region.GroupName))
                                {
                                    continue;
                                }
                            }
                            //it's legal to have Start > End for the case (mitochondria) where an interval spans the "end" of a circular chromosome
                            if ((region.Start < region.End) && region.Start + region.StartProbeLength + region.EndProbeLength > region.End)
                            {
                                throw new Exception(string.Format(
                                    "Error: Region {0} in manifest {1} is not valid: Probe lengths are larger than the start...end interval.",
                                    region.Name, manifestName));
                            }
                            Regions.Add(region);
                        }
                    }
                }
            }
            if (lineNumber == -1)
                throw new Exception($"Invalid manifest: No [Regions] section seen in {manifestName}");
            if (Regions.Count == 0)
                throw new Exception($"Error: No regions seen in manifest {manifestName}");

            // Sort regions:
            Regions.Sort();

            // If ExcludeRegionGroups list non-existent region groups, that's potentially a problem, since the user may have 
            // made a typo and attempted to exclude a group.  Log a warning but continue:
            if (ExcludeRegionGroups != null)
            {
                foreach (string ExcludeGroup in ExcludeRegionGroups)
                {
                    if (!GroupNamesSeen.Contains(ExcludeGroup))
                    {
                        if (Error != null)
                        {
                            Error(
                                $"Warning: ExcludeRegionGroups setting specified group '{ExcludeGroup}', " +
                                "but no regions in manifest file {manifestName} are assigned to this group;" +
                                "no region exclusion performed");
                        }
                    }
                }
            }

        }


        /// <summary>
        /// Get the sorted manifest regions for each chromosome.
        /// </summary>
        /// <returns>Dictionary of chromosome name to sorted manifest regions of the chromosome.</returns>
        public Dictionary<string, List<ManifestRegion>> GetManifestRegionsByChromosome()
        {
            Dictionary<string, List<ManifestRegion>> regionsByChrom = new Dictionary<string, List<ManifestRegion>>();

            foreach (var region in Regions)
            {
                if (!regionsByChrom.ContainsKey(region.Chromosome)) { regionsByChrom[region.Chromosome] = new List<ManifestRegion>(); }
                regionsByChrom[region.Chromosome].Add(region);
            }

            foreach (string chrom in regionsByChrom.Keys)
            {
                regionsByChrom[chrom].Sort();
            }

            return regionsByChrom;
        }

        public class ManifestRegion : IComparable<ManifestRegion>
        {
            public string Chromosome;
            public int End;
            public int EndProbeLength = 0;
            public string Name;
            public int Start;
            public int StartProbeLength = 0;
            public string GroupName; // Regions for the same gene (or other locus) can be grouped together, then excluded by group. 

            public ManifestRegion() { }

            //for cloning
            public ManifestRegion(ManifestRegion region)
            {
                Chromosome = region.Chromosome;
                End = region.End;
                EndProbeLength = region.EndProbeLength;
                Name = region.Name;
                Start = region.Start;
                StartProbeLength = region.StartProbeLength;
                GroupName = region.GroupName;
            }

            /// <summary>
            /// Sort by chromosome name (asciibetical), then start position, then end position.
            /// If region[n] starts after a variant, so will regions with index >n
            /// But note that if region[n] ends before a variant, regions <n may still cover the variant
            /// </summary>
            public int CompareTo(ManifestRegion other)
            {
                if (Chromosome != other.Chromosome) return string.Compare(Chromosome, other.Chromosome);
                if (Start < other.Start) return -1;
                if (Start > other.Start) return 1;
                if (End < other.End) return -1;
                if (End > other.End) return 1;
                return string.Compare(Name, other.Name);
            }

            /// <summary>
            /// Sort by chromosome index from fasta file, then start position, then end position.
            /// If region[n] starts after a variant, so will regions with index >n
            /// But note that if region[n] ends before a variant, regions <n may still cover the variant
            /// </summary>
            public int CompareTo(ManifestRegion other, Dictionary<string, int> chromosomeIndexLookup)
            {
                int compareChromosome = chromosomeIndexLookup[Chromosome].CompareTo(chromosomeIndexLookup[other.Chromosome]);
                if (compareChromosome != 0) return compareChromosome;
                if (Start < other.Start) return -1;
                if (Start > other.Start) return 1;
                if (End < other.End) return -1;
                if (End > other.End) return 1;
                return string.Compare(Name, other.Name);
            }

            public TargetInterval GetTargetInterval(int refIndex = -1)
            {
                // remove the probe portions from the target interval
                int start = Start + StartProbeLength;
                int end = End - EndProbeLength;

                return new TargetInterval(Chromosome, refIndex, start, end);
            }
        }
    }

    // this class defines the regions specified in the manifest file
    public class TargetInterval : IComparable<TargetInterval>
    {
        #region members

        public int Begin;
        public int End;
        private bool _hasHashCode;
        private int _hashCode;
        public int ReferenceIndex;
        public string ReferenceName;

        #endregion

        // constructor
        public TargetInterval(string referenceName, int referenceIndex, int begin, int end)
        {
            ReferenceName = referenceName;
            ReferenceIndex = referenceIndex;
            Begin = begin;
            End = end;
        }

        public int CompareTo(TargetInterval ti)
        {
            // If other is not a valid object reference, this instance is greater.
            if (ti == null) return 1;

            // The Region comparison depends on the comparison of 
            // the reference name, begin, and end
            if (ReferenceIndex == ti.ReferenceIndex)
            {
                if (Begin == ti.Begin) return End.CompareTo(ti.End);
                return Begin.CompareTo(ti.Begin);
            }

            return ReferenceIndex.CompareTo(ti.ReferenceIndex);
        }

        public override int GetHashCode()
        {
            if (!_hasHashCode)
            {
                _hashCode = Begin.GetHashCode() ^ End.GetHashCode() ^ ReferenceIndex.GetHashCode();
                _hasHashCode = true;
            }

            return _hashCode;
        }

        public override bool Equals(object obj)
        {
            if (obj == null) return false;
            TargetInterval other = (TargetInterval)obj;
            if (ReferenceIndex != other.ReferenceIndex) return false;
            if (Begin != other.Begin) return false;
            if (End != other.End) return false;
            return true;
        }

        /// <summary>
        ///     returns true if this region overlaps with the specified region
        /// </summary>
        public bool Overlaps(TargetInterval ti)
        {
            if (ReferenceIndex != ti.ReferenceIndex) return false;
            if ((End < ti.Begin) || (Begin > ti.End)) return false;
            return true;
        }
    }

    public class IntervalSet
    {
        #region members
        public int GenomeIndex;
        public string ManifestName;
        public string ReferencePath;
        public HashSet<TargetInterval> TargetIntervals;
        #endregion

        // constructor
        public IntervalSet()
        {
            TargetIntervals = new HashSet<TargetInterval>();
        }
    }
}
