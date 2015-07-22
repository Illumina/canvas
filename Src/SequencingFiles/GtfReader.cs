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

        public class GtfEntry: IComparable<GtfEntry>
        {
            public string seqname;
            public string source;
            public string feature;
            public int start;
            public int end;
            public float score;
            public Strand strand;
            public int frame;
            public string geneID;
            public Dictionary<string, string> attributes;

            // logic taken from ManifestRegion.CompareTo
            public int CompareTo(GtfEntry other) 
            {
                if (seqname != other.seqname) return seqname.CompareTo(other.seqname);
                if (start < other.start) return -1;
                if (start > other.start) return 1;
                if (end < other.end) return -1;
                if (end > other.end) return 1;
                return geneID.CompareTo(other.geneID);
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
            if (requiredAttributes == null) { requiredAttributes = new List<string>(); }
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

                    // Parse frame
                    // May also be something other than a number. Default to 0.
                    int frame;
                    if (!int.TryParse(lineSplit[7], out frame)) { frame = 0; }
                    
                    // Parse attributes
                    var attribs = ParseAttributes(lineSplit[8]);
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

                    // Construct the GtfEntry
                    GtfEntry entry = new GtfEntry()
                    {
                        seqname = seqname,
                        source = source,
                        feature = feature,
                        start = start,
                        end = end,
                        score = score,
                        strand = strand,
                        geneID = attribs.ContainsKey("gene_id") ? attribs["gene_id"] : null,
                    };
                    if (storeAttributes) entry.attributes = attribs;
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
            char sep;
            if (attributes.Contains('=')) 
            {
                sep = '=';
            }
            else
            {
                sep = ' ';
            }

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
            foreach (GtfEntry entry in gtfEntries)
            {
                if (entry.feature.ToLower() == "exon") {
                    yield return entry;
                }
            }
        }

        public List<Transcript> GetTranscripts(IEnumerable<GtfEntry> gtfEntries)
        {
            Dictionary<string, Transcript> transcriptByID = new Dictionary<string, Transcript>();
            foreach (GtfEntry entry in gtfEntries) 
            {
                if (!entry.attributes.ContainsKey("transcript_id")) { continue; }
                string id = entry.attributes["transcript_id"];

                Transcript transcript;
                if (!transcriptByID.ContainsKey(id))
                {
                    transcript = new Transcript()
                    {
                        chrom = entry.seqname,
                        start = entry.start,
                        stop = entry.end,
                        strand = entry.strand,
                        id = id,
                        geneID = entry.geneID
                    };
                    transcriptByID[id] = transcript;
                }
                else 
                {
                    transcript = transcriptByID[id];
                    if (transcript.start > entry.start) { transcript.start = entry.start; }
                    if (transcript.stop < entry.end) { transcript.stop = entry.end; }
                }

                switch (entry.feature.ToLower()) 
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
}
