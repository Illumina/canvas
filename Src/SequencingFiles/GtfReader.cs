using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace SequencingFiles
{
    /// <summary>
    /// Parser for GTF annotation files
    /// </summary>
    public class GtfReader
    {
        public enum Strand { plus, minus, unknown };

        public static Strand ReverseStrand(Strand strand)
        {
            switch (strand)
            {
                case Strand.plus: return Strand.minus;
                case Strand.minus: return Strand.plus;
                case Strand.unknown: return Strand.unknown;
                default: return Strand.unknown;
            }
        }

        /// <summary>
        /// Represents one entry (line) from the GTF file, e.g. one exon
        /// </summary>
        public class GtfEntry: IComparable<GtfEntry>
        {
            // From mandatory GTF fields
            public readonly string Seqname;
            public readonly string Source;
            public readonly string Feature;
            public readonly int Start;
            public readonly int End;
            public readonly float Score;
            public readonly Strand Strand;
            public readonly int? Frame; // null if not a number
            // From extra attributes; maybe null 
            public string GeneID; // should be unique
            private string geneName; // clear text name; may not be unique
            private Dictionary<string, string> attributes;

            public GtfEntry(string seqname, string source, string feature, int start, int end,
                float score, Strand strand, int? frame, Dictionary<string, string> attributes=null, bool storeAttributes=true)
            {
                this.Seqname = seqname;
                this.Source = source;
                this.Feature = feature;
                this.Start = start;
                this.End = end;
                this.Score = score;
                this.Strand = strand;
                this.Frame = frame;  
                if (attributes != null)
                {
                    attributes.TryGetValue("gene_id", out GeneID);
                    attributes.TryGetValue("gene_name", out geneName);
                }
                if (storeAttributes) this.attributes = attributes;
            }

            // logic taken from ManifestRegion.CompareTo
            public int CompareTo(GtfEntry other) 
            {
                if (Seqname != other.Seqname) return Seqname.CompareTo(other.Seqname);
                if (Start < other.Start) return -1;
                if (Start > other.Start) return 1;
                if (End < other.End) return -1;
                if (End > other.End) return 1;
                return GeneID.CompareTo(other.GeneID);
            }

            /// <summary>
            /// Return the value of one extra GTF attribute for this entry
            /// </summary>
            /// <param name="attribute">Name of the attribute to lookup</param>
            /// <returns>Attribute value as a string or null if not present (or attributes) are not stored</returns>
            public string GetAttribute(string attribute)
            {
                if ((attributes != null) && attributes.ContainsKey(attribute))
                    return attributes[attribute];
                return null;
            }
            /// <summary>
            /// Check if this GTF entry contains a given attribute
            /// </summary>
            /// <param name="attribute">Name of the attribute to lookup</param>
            /// <returns>True if the attribute is present, otherwise false</returns>
            public bool HasAttribute(string attribute)
            {
                return ((attributes != null) && attributes.ContainsKey(attribute));
            }
            /// <summary>
            /// Return the gene name of this element, based on the "gene_name" attribute if present or "gene_id" otherwise.
            /// Gene_name may not be unique.
            /// Null if neither is present.
            /// </summary>
            public string GeneName
            { get
                { return geneName ?? GeneID; }
            }
        }

        public readonly string gtfFile;
        
        // For logging, we need a lambda that accepts a string and returns nothing (e.g. Console.WriteLine(string))
        // The declared type of such a function uses the built-in delegate Action<string>
        Action<string> errorHandler, logHandler;

        public List<string> comments { get; private set; }

        public class Transcript
        {
            public string chrom;
            public int start, stop; // GTF coordinates are 1-based
            public Strand strand;
            public string id;
            public string geneID;
            public List<GtfEntry> exons = new List<GtfEntry>();
            public GtfEntry startCodon = null;
            public GtfEntry stopCodon = null;
            public List<GtfEntry> CDSs = new List<GtfEntry>();
            public List<GtfEntry> otherFeatures = new List<GtfEntry>();
            public int CDSStart 
            {
                get 
                {
                    if (CDSs.Any())
                    {
                        return CDSs.First().Start;
                    }
                    else
                    {
                        return stop + 1; // used by the refFlat format
                    }
                }
            }
            public int CDSStop
            {
                get 
                {
                    if (CDSs.Any())
                    {
                        return CDSs.Last().End;
                    }
                    else 
                    {
                        return stop; // used by the refFlat format
                    }
                }
            }
            public string RefFlat
            {
                /// GTF: 1-based coordinates
                /// refFlat: 0-based coordinates
                get 
                {
                    string exonStarts = String.Join(",", exons.Select(e => (e.Start - 1).ToString())) + ",";
                    string exonEnds = String.Join(",", exons.Select(e => e.End.ToString())) + ",";
                    return String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}",
                        geneID, id, chrom, strand.ToString2(), start - 1, stop,
                        CDSStart - 1, CDSStop, exons.Count, exonStarts, exonEnds);
                }
            }
        }
        
        public GtfReader(string gtfFilename, 
            Action<string> log = null, Action<string> error = null)
        {
            errorHandler = error == null ? Console.Error.WriteLine : error;
            logHandler = log == null? Console.WriteLine: log;
            gtfFile = gtfFilename;
            comments = new List<string>();
        }

        // Returns an iterator over all GTF entries in the file
        public IEnumerable<GtfEntry> parseGtfFile(bool storeAttributes,
            List<string> requiredAttributes = null, bool ignoreUnknownStrand = true)
        {
            comments.Clear();
            // TODO - handle "track" entries.
            using (StreamReader reader = new StreamReader(gtfFile))
            {
                while (!reader.EndOfStream)
                {
                    string fileLine = reader.ReadLine();
                    if (fileLine.StartsWith("#")) 
                    {
                        comments.Add(fileLine);
                        continue; // header comments
                    }

                    string[] lineSplit = fileLine.Split('\t');
                    if (lineSplit.Length != 9)
                        throw new ApplicationException("Invalid GTF line:\n" + fileLine);

                    // Parse seqname
                    string seqname = lineSplit[0];

                    // Parse source
                    string source = lineSplit[1];

                    // Parse feature
                    string feature = lineSplit[2];

                    // Parse start, end
                    int start = int.Parse(lineSplit[3]);
                    int end = int.Parse(lineSplit[4]);
                    if (start > end)
                    {
                        errorHandler("Start > Stop: " + fileLine);
                        continue;
                    }

                    // Parse score
                    // May have something other than a number - in which case, use 0.
                    float score;
                    if (!float.TryParse(lineSplit[5], out score)) { score = 0; }
                    
                    // Parse strand
                    var strand = ParseStrand(lineSplit[6]);
                    if (strand == Strand.unknown && ignoreUnknownStrand)
                    {
                        errorHandler("Unknown strand: " + fileLine);
                        continue;
                    }

                    // Parse frame. 'null' if something other than a number.
                    int tmpFrame;
                    int? frame = null;
                    if (int.TryParse(lineSplit[7], out tmpFrame))
                        frame = tmpFrame;
                    
                    // Parse attributes
                    var attribs = ParseAttributes(lineSplit[8]);
                    if (requiredAttributes != null)
                    {
                        bool hasAllRequiredAttribs = true;
                        List<string> missingAttribs = new List<string>();
                        foreach (string attrib in requiredAttributes)
                        {
                            if (!attribs.ContainsKey(attrib))
                            {
                                hasAllRequiredAttribs = false;
                                missingAttribs.Add(attrib);
                            }
                        }
                        if (!hasAllRequiredAttribs)
                        {
                            errorHandler(String.Format("Missing {0}: {1}", String.Join(", ", missingAttribs), fileLine));
                            continue;
                        }
                    }
                    // Construct the GtfEntry
                    GtfEntry entry = new GtfEntry(seqname, source, feature, start, end, score, strand, frame,
                        attribs, storeAttributes);
                    yield return entry;
                }
            }
        }

        private Strand ParseStrand(string strand)
        {
            switch (strand)
            {
                case "+":
                    return Strand.plus;
                case "-":
                    return Strand.minus;
                default:
                    return Strand.unknown;
            }
        }
        
        // Parse attributes of an GTF entry
        private Dictionary<string, string> ParseAttributes(string attributes)
        {
            // Determine the key/value separator character. May be ' ' or '=' 
            // If '=' is present, we'll use it. Otherwise we'll use ' '.
            char sep = (attributes.Contains('=') ? '=' : ' ');
            var attrs = new Dictionary<string, string>();
            foreach (var item in attributes.Split(';')) // Attrs. seperated by '; '
            {
                if (item == "") continue;
                var attr = item.Trim().Split(sep); // Split Key / Value
                string key, value;
                if (attr.Length != 2)
                {
                    key = attr[0];
                    value = "";
                }
                else
                {
                    key = attr[0];
                    value = attr[1].Trim('"'); // Trim '"' around value - present in some formats.
                }
                attrs[key] = value; 
            }
            return attrs;
        }

        public IEnumerable<GtfEntry> GetExons(IEnumerable<GtfEntry> gtfEntries)
        {
            return gtfEntries.Where(e => e.Feature.ToLower() == "exon");
        }

        public List<Transcript> GetTranscripts(IEnumerable<GtfEntry> gtfEntries)
        {
            Dictionary<string, Transcript> transcriptByID = new Dictionary<string, Transcript>();
            foreach (GtfEntry entry in gtfEntries) 
            {
                string id = entry.GetAttribute("transcript_id");
                if (id == null) continue;
                Transcript transcript;
                if (!transcriptByID.ContainsKey(id))
                {
                    transcript = new Transcript()
                    {
                        chrom = entry.Seqname,
                        start = entry.Start,
                        stop = entry.End,
                        strand = entry.Strand,
                        id = id,
                        geneID = entry.GeneID
                    };
                    transcriptByID[id] = transcript;
                }
                else 
                {
                    transcript = transcriptByID[id];
                    if (transcript.start > entry.Start) { transcript.start = entry.Start; }
                    if (transcript.stop < entry.End) { transcript.stop = entry.End; }
                }

                switch (entry.Feature.ToLower()) 
                {
                    case "exon":
                        transcript.exons.Add(entry);
                        break;
                    case "start_codon":
                        transcript.startCodon = entry;
                        break;
                    case "stop_codon":
                        transcript.stopCodon = entry;
                        break;
                    case "cds":
                        transcript.CDSs.Add(entry);
                        break;
                    default:
                        transcript.otherFeatures.Add(entry);
                        break;
                }
            }

            foreach (Transcript transcript in transcriptByID.Values)
            {
                transcript.exons.Sort();
                transcript.CDSs.Sort();
                transcript.otherFeatures.Sort();
            }

            return transcriptByID.Values.ToList();
        }
    }

    public static class StrandExtensions
    {
        public static string ToString2(this GtfReader.Strand s)
        {
            switch (s)
            {
                case GtfReader.Strand.plus:
                    return "+";
                case GtfReader.Strand.minus:
                    return "-";
                case GtfReader.Strand.unknown:
                default:
                    return ".";
            }
        }
    }
}
