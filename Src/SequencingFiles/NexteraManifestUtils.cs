using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using SequencingFiles;

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
    }
}