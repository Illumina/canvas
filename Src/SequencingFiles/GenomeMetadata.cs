using System;
using System.Collections.Generic;
using System.IO;
using System.Security.Cryptography;
using System.Text;
using System.Text.RegularExpressions;
using System.Xml;
using ProtoBuf;

namespace SequencingFiles
{
    public class GenomeMetadata
    {
        #region Members
        public const int BufferSize = 10485760; // 1024*1024*10
        protected static Regex IndexRegex;
        public long Length { get; set; }
        public string Name { get; set; }
        public List<SequenceMetadata> Sequences { get; set; }
        public long KnownBases { get; set; }
        public string Species
        {
            get
            {
                if (Sequences.Count == 0)
                    return null;
                return Sequences[0].Species;
            }
        }
        public string Build
        {
            get
            {
                if (Sequences.Count == 0)
                    return null;
                return Sequences[0].Build;
            }
        }
        #endregion

        public List<SequenceMetadata> GetChromosomesIncludingNull()
        {
            List<SequenceMetadata> list = new List<SequenceMetadata>();
            list.AddRange(Sequences);
            list.Add(null);
            return list;
        }

        // constructor
        public enum GenomeFolderState
        {
            Ready = 0, // No action needed
            RequireWritableFolder, // Ensure we have a writable genome folder, but do nothing else (yet)
            RequireImport, // Must import (in a writable folder)
            RequireFASTACombine, // Must combine FASTA files (in a writable folder) and import
        }

        public GenomeMetadata()
        {
            Sequences = new List<SequenceMetadata>();
            IndexRegex = new Regex(@"^(\d+)\t>(.+)$", RegexOptions.Compiled);
        }

        /// <summary>
        /// Adds a new reference sequence to the Sequences list
        /// </summary>
        private void AddReferenceSequence(string name, long length, string fastaPath, FastaHeaderMetadata header,
            ref HashSet<string> referenceNames, long knownBaseLength)
        {
            // sanity check
            if (referenceNames.Contains(name))
            {
                throw new ApplicationException(string.Format(
                    "An attempt was made to load two reference sequences with the same exact names ({0})", name));
            }
            referenceNames.Add(name);

            // add the reference sequence
            SequenceMetadata sequence = new SequenceMetadata
            {
                FastaPath = fastaPath,
                Length = length,
                Name = name,
                Build = header.Build,
                Checksum = header.Checksum,
                IsCircular = header.IsCircular,
                Ploidy = header.Ploidy,
                Species = header.Species,
                KnownBases = knownBaseLength
            };

            // update species and build from fasta path if in iGenomes format
            IGenomesReferencePath iGenomesReference = IGenomesReferencePath.GetReferenceFromFastaPath(sequence.FastaPath);
            if (iGenomesReference != null)
            {
                if (string.IsNullOrEmpty(sequence.Build))
                    sequence.Build = iGenomesReference.Build;
                if (string.IsNullOrEmpty(sequence.Species))
                    sequence.Species = iGenomesReference.Species;
            }

            Sequences.Add(sequence);
        }

        static public SequenceType ParseSequenceType(string type)
        {
            if (type == null) return SequenceType.Unknown;
            switch (type.ToLowerInvariant())
            {
                case "althaplotype":
                    return SequenceType.AltHaplotype;
                case "autosome":
                    return SequenceType.Autosome;
                case "contig":
                    return SequenceType.Contig;
                case "decoy":
                    return SequenceType.Decoy;
                case "mitochondria":
                    return SequenceType.Mitochondria;
                case "sex":
                case "allosome":
                    return SequenceType.Allosome;
                case "other":
                    return SequenceType.Other;
                default:
                    return SequenceType.Unknown;
            }
        }


        /// <summary>
        /// Populates the genome metadata from an XML file
        /// </summary>
        public void Deserialize(string inputFilename)
        {
            // open the XML file
            inputFilename = Path.GetFullPath(inputFilename);
            string directory = Path.GetDirectoryName(inputFilename);
            Length = 0;
            KnownBases = 0; // initial
            int refIndex = 0;
            IGenomesReferencePath iGenomesReference = IGenomesReferencePath.GetReferenceFromFastaPath(directory);

            // use StreamReader to avoid URI parsing of filename that will cause problems with 
            // certain characters in the path (#). 
            using (var streamReader = new StreamReader(inputFilename))
            using (var xmlReader = XmlReader.Create(streamReader))
            {
                while (xmlReader.Read())
                {
                    XmlNodeType nType = xmlReader.NodeType;

                    // handle 
                    if (nType == XmlNodeType.Element)
                    {
                        // retrieve the genome variables
                        if (xmlReader.Name == "sequenceSizes")
                        {
                            Name = xmlReader.GetAttribute("genomeName");
                            if (iGenomesReference != null && string.IsNullOrEmpty(Name))
                                Name = iGenomesReference.ToString();
                        }

                        // retrieve the chromosome variables
                        if (xmlReader.Name == "chromosome")
                        {
                            SequenceMetadata refSeq = new SequenceMetadata
                            {
                                FastaPath = Path.Combine(directory, xmlReader.GetAttribute("fileName")),
                                Name = xmlReader.GetAttribute("contigName"),
                                Index = refIndex++,
                                Length = long.Parse(xmlReader.GetAttribute("totalBases")),
                                Type = ParseSequenceType(xmlReader.GetAttribute("type"))
                            };
                            Length += refSeq.Length;

                            refSeq.Build = xmlReader.GetAttribute("build");
                            refSeq.Species = xmlReader.GetAttribute("species");

                            // update species and build from fasta path if in iGenomes format
                            if (iGenomesReference != null)
                            {
                                if (string.IsNullOrEmpty(refSeq.Build))
                                    refSeq.Build = iGenomesReference.Build;
                                if (string.IsNullOrEmpty(refSeq.Species))
                                    refSeq.Species = iGenomesReference.Species;
                            }

                            string isCircular = xmlReader.GetAttribute("isCircular");
                            if (!string.IsNullOrEmpty(isCircular))
                                refSeq.IsCircular = (isCircular == "true");

                            string ploidy = xmlReader.GetAttribute("ploidy");
                            if (!string.IsNullOrEmpty(ploidy)) refSeq.Ploidy = int.Parse(ploidy);

                            string md5 = xmlReader.GetAttribute("md5");
                            if (!string.IsNullOrEmpty(md5)) refSeq.Checksum = md5;

                            string knownBases = xmlReader.GetAttribute("knownBases");
                            if (!string.IsNullOrEmpty(knownBases))
                            {
                                refSeq.KnownBases = long.Parse(knownBases);
                                KnownBases += refSeq.KnownBases;
                            }

                            Sequences.Add(refSeq);
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Scans the reference sequence list and returns the specified sequence metadata
        /// TODO: nulls are evil. throw an exception rather than return null. have them use TryGetSequence instead
        /// TODO: create lookup table to make this faster? 
        /// </summary>
        public SequenceMetadata GetSequence(string sequenceName)
        {
            foreach (SequenceMetadata sequence in Sequences)
            {
                if (string.Equals(sequence.Name, sequenceName, StringComparison.OrdinalIgnoreCase))
                    return sequence;
            }
            return null;
        }

        /// <summary>
        /// Scans the reference sequence list and returns the specified sequence metadata if found
        /// TODO: create lookup table to make this faster? 
        /// </summary>
        public bool TryGetSequence(string sequenceName, out SequenceMetadata foundSequence)
        {
            foreach (SequenceMetadata sequence in Sequences)
            {
                if (string.Equals(sequence.Name, sequenceName, StringComparison.OrdinalIgnoreCase))
                {
                    foundSequence = sequence;
                    return true;
                }
            }

            foundSequence = null;
            return false;
        }

        /// <summary>
        /// Serializes the genome metadata to an XML file
        /// </summary>
        public void Serialize(string outputFilename)
        {
            // open the XML file
            // Initialize with StreamWriter to avoid URI parsing of input filename which doesn't 
            // like # in the path
            using (XmlTextWriter xmlWriter = new XmlTextWriter(new StreamWriter(outputFilename, false, Encoding.ASCII)))
            {
                xmlWriter.Formatting = Formatting.Indented;
                xmlWriter.IndentChar = '\t';
                xmlWriter.Indentation = 1;

                // write all of our sequences
                xmlWriter.WriteStartElement("sequenceSizes");
                if (!string.IsNullOrEmpty(Name)) xmlWriter.WriteAttributeString("genomeName", Name);

                foreach (SequenceMetadata refSeq in Sequences)
                {
                    xmlWriter.WriteStartElement("chromosome");

                    // required for compatibility with CASAVA 1.8
                    xmlWriter.WriteAttributeString("fileName", Path.GetFileName(refSeq.FastaPath));
                    xmlWriter.WriteAttributeString("contigName", refSeq.Name);
                    xmlWriter.WriteAttributeString("totalBases", refSeq.Length.ToString());

                    // additional attributes for MiSeq
                    if (!string.IsNullOrEmpty(refSeq.Build)) xmlWriter.WriteAttributeString("build", refSeq.Build);
                    xmlWriter.WriteAttributeString("isCircular", (refSeq.IsCircular ? "true" : "false"));
                    if (!string.IsNullOrEmpty(refSeq.Checksum)) xmlWriter.WriteAttributeString("md5", refSeq.Checksum);
                    xmlWriter.WriteAttributeString("ploidy", refSeq.Ploidy.ToString());
                    if (!string.IsNullOrEmpty(refSeq.Species))
                        xmlWriter.WriteAttributeString("species", refSeq.Species);
                    xmlWriter.WriteAttributeString("knownBases", refSeq.KnownBases.ToString());
                    xmlWriter.WriteAttributeString("type", refSeq.Type.ToString());

                    xmlWriter.WriteEndElement();
                }
                xmlWriter.WriteEndElement();

                // close the XML file
                xmlWriter.Close();
            }
        }

        /// <summary>
        ///     Retrieves the FASTA filenames from the specified directory
        /// </summary>
        /// <returns>A list of FASTA filenames</returns>
        private static List<string> GetFastaFilenames(string directory)
        {
            List<string> fastaFilenames = new List<string>();

            if (Directory.Exists(directory))
            {
                DirectoryInfo info = new DirectoryInfo(directory);
                foreach (FileInfo fi in info.GetFiles("*.fa")) fastaFilenames.Add(fi.FullName);
                foreach (FileInfo fi in info.GetFiles("*.fasta")) fastaFilenames.Add(fi.FullName);
            }

            return fastaFilenames;
        }

        /// <summary>
        /// Checks if a file exists and was written to after a reference timepoint (FASTA write time)
        /// </summary>
        /// <returns>true if the file exists and the file is more recent than the FASTA file</returns>
        private static bool CheckFile(string filePath, DateTime compareTime)
        {
            bool isFileGood = File.Exists(filePath);
            if (isFileGood)
            {
                DateTime fileWriteTime = File.GetLastWriteTime(filePath);
                if (fileWriteTime < compareTime)
                {
                    // Update 6/3/13: Don't be strict about the modification times.  In practice, iGenomes can't guarantee (especially
                    // after references are copied here and there) that the modification times will be preserved.  We're more likely
                    // to create problems by rejecting index files with older modification times than we are to prevent issues due to
                    // editing the FASTA file.
                    Console.WriteLine("Warning: Modification time of FASTA file is more recent than {0}.  If FASTA file contents have been modified, please re-generate indexes to ensure they are valid.", filePath);
                }
            }
            return isFileGood;
        }

        /// <summary>
        ///     Determines the state of a reference genome folder - is it ready to go, do we need to double-check that the
        ///     folder is writable, do we need to import, do we need to combine FASTA files.
        /// </summary>
        public static GenomeFolderState CheckReferenceGenomeFolderState(string directory, bool requireBWT, bool requireBowtie)
        {
            List<string> fastaFilenames = GetFastaFilenames(directory);
            DateTime mostRecentFastaFile = DateTime.MinValue;

            // If ther's more than one FASTA, then we need to import.
            if (fastaFilenames.Count > 1)
            {
                Console.WriteLine(">>> Multiple FASTA files -> require import!");
                return GenomeFolderState.RequireImport;
            }
            if (fastaFilenames.Count == 0)
            {
                throw new Exception(string.Format("Error: No reference genome FASTA file (genome.fa) found in folder {0}", directory));
            }

            bool requireWritableFolder = false;

            // check the derivative files
            foreach (string fastaPath in fastaFilenames)
            {
                // retrieve the FASTA time
                DateTime fastaWriteTime = File.GetLastWriteTime(fastaPath);
                if (fastaWriteTime > mostRecentFastaFile) mostRecentFastaFile = fastaWriteTime;

                if (requireBowtie)
                {
                    foreach (string suffix in new[] { ".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt" })
                    {
                        if (!CheckFile(Path.Combine(Path.GetDirectoryName(fastaPath),
                            Path.GetFileNameWithoutExtension(fastaPath)) + suffix, fastaWriteTime))
                        {
                            Console.WriteLine(">>> Require bowtie -> require writable");
                            requireWritableFolder = true;
                            break;
                        }
                    }
                }
                if (requireBWT)
                {
                    // Check the canonical iGenomes path:
                    string tempFolder = Path.GetDirectoryName(directory);
                    tempFolder = Path.Combine(tempFolder, "BWAIndex");
                    string fastaFileName = Path.GetFileName(fastaPath);
                    bool HasBWT = CheckFile(Path.Combine(tempFolder, fastaFileName + ".bwt"), fastaWriteTime);
                    if (!HasBWT)
                    {
                        // Check in the FASTA folder too:
                        HasBWT = CheckFile(fastaPath + ".bwt", fastaWriteTime);
                    }
                    if (!HasBWT)
                    {
                        Console.WriteLine(">>> Require bwt -> require writable");
                        requireWritableFolder = true;
                    }
                }
                // Check for FOO.fa.fai and FOO.dict:
                if (!CheckFile(fastaPath + ".fai", fastaWriteTime))
                {
                    Console.WriteLine(">>> FAI file missing -> import needed");
                    return GenomeFolderState.RequireImport;
                }
                string dictFilename = Path.Combine(Path.GetDirectoryName(fastaPath), Path.GetFileNameWithoutExtension(fastaPath)) + ".dict";
                if (!CheckFile(dictFilename, fastaWriteTime))
                {
                    Console.WriteLine(">>> Dict file missing -> import needed");
                    return GenomeFolderState.RequireImport;
                }
            }

            // check for genomesize.xml file
            if (!CheckFile(Path.Combine(directory, "GenomeSize.xml"), mostRecentFastaFile))
            {
                Console.WriteLine(">>> GenomeSize.xml file missing -> import needed");
                return GenomeFolderState.RequireImport;
            }

            if (requireWritableFolder) return GenomeFolderState.RequireWritableFolder;
            return GenomeFolderState.Ready;
        }

        public void ImportFromFastaFiles(string directory)
        {
            // initialize
            directory = Path.GetFullPath(directory);
            List<string> fastaFilenames = GetFastaFilenames(directory);
            HashSet<string> referenceNames = new HashSet<string>();

            // Index the FASTA file:
            foreach (string fastaPath in fastaFilenames)
            {
                List<FastaHeaderMetadata> metadataList = new List<FastaHeaderMetadata>();

                List<RefIndexEntry> referenceIndexes = new List<RefIndexEntry>();
                ReferenceIndexer.CreateFASTAIndexFiles(fastaPath, metadataList, referenceIndexes);
                for (int chromosomeIndex = 0; chromosomeIndex < metadataList.Count; chromosomeIndex++)
                {
                    AddReferenceSequence(referenceIndexes[chromosomeIndex].Name, referenceIndexes[chromosomeIndex].Length,
                        Path.Combine(directory, fastaPath), metadataList[chromosomeIndex], ref referenceNames,
                        referenceIndexes[chromosomeIndex].KnownBasesLength);
                }
            }

            // set the genome information
            Name = Path.GetFileName(directory);
            if (Name.ToLowerInvariant() == "chromosomes" || Name.ToLowerInvariant() == "wholegenomefasta")
            {
                // Try to give a unique name rather than generic "chromosomes" or "wholegenomefasta":
                Name = Path.GetDirectoryName(directory);
                int position = Name.ToLowerInvariant().IndexOf("igenomes");
                if (position != -1)
                {
                    Name = Name.Substring(position + 8);
                    if (Name[0] == '\\') Name = Name.Substring(1);
                }
            }

            Length = 0;
            KnownBases = 0;
            foreach (SequenceMetadata refSeq in Sequences)
            {
                Length += refSeq.Length;
                KnownBases += refSeq.KnownBases;
            }
        }

        public enum SequenceType : int
        {
            Unknown = 0, // not classified. Default for older GenomeSize.xml files where sequence type is not specified
            Autosome, // main chromsomes (chr1, chr2...)
            Mitochondria, // chrM
            Allosome, // sex chromosome (chrX, chrY)
            Contig, // Unplaced or unlocalized contigs (e.g. chr1_KI270711v1_random, chrUn_KI270530v1)
            Decoy,
            AltHaplotype,
            Other, // currently only chrEBV has this classification
        }

        [Serializable]
        [ProtoContract]
        public class SequenceMetadata : IComparable<SequenceMetadata>
        {
            #region member variables
            [ProtoMember(1)]
            public bool IsCircular = false;
            [ProtoMember(2)]
            public int Index { get; set; }
            [ProtoMember(3)]
            public int Ploidy { get; set; }
            [ProtoMember(4)]
            public long Length { get; set; }
            [ProtoMember(5)]
            public string Build { get; set; }
            [ProtoMember(6)]
            public string Checksum { get; set; }
            [ProtoMember(7)]
            public string FastaPath { get; set; }
            [ProtoMember(8)]
            public string Name { get; set; }
            [ProtoMember(9)]
            public string Species { get; set; }
            [ProtoMember(10)]
            public long KnownBases { get; set; } // bases that are not 'N'
            [ProtoMember(11)]
            public SequenceType Type { get; set; }

            #endregion

            #region constructors

            public SequenceMetadata()
            {
                Ploidy = 2;
                Index = -1;
            }

            public SequenceMetadata(string name, long referenceLength, int index)
                : this()
            {
                Name = name;
                Length = referenceLength;
                Index = index;
            }

            #endregion

            /// <summary>
            /// 
            /// TODO: completely eliminate the use of this method. All human GenomeSize.xml will have the type specified
            /// 
            /// Checks if the chromosome is an autosome.
            /// Assumes species is Homo_sapiens 
            /// This metadata should be part of the GenomeSize.xml and come from the Type member, but we fall back
            /// to checking the string if the attribute is missing.
            /// </summary>
            /// <param name="chr">Name of chromosome to check.</param>
            /// <returns></returns>
            public static bool IsAutosome(string chr)
            {
                string id = chr.Replace("chr", "");
                double number;
                bool isNumber = double.TryParse(id, out number);
                return isNumber;
            }

            public bool IsAutosome()
            {
                switch (this.Type)
                {
                    case SequenceType.Autosome:
                        return true;
                    case SequenceType.Contig:
                    case SequenceType.Decoy:
                    case SequenceType.Mitochondria:
                    case SequenceType.Allosome:
                    case SequenceType.Other:
                        return false;
                    case SequenceType.Unknown:
                    default:
                        return IsAutosome(Name);
                }
            }

            /// <summary>
            /// TODO: DO NOT use this. There is too much ambiguity here when chrom is not found. The caller should use GenomeMetadata.TryGetSequence instead and then call IsDecoyOrOther on the returned SequenceMetadata. Don't have the GenomeMetadata? Tell the caller to give it to you.
            /// </summary>
            /// <param name="chrom"></param>
            /// <param name="chromosomes"></param>
            /// <returns></returns>
            public static bool IsDecoyOrOther(string chrom, List<GenomeMetadata.SequenceMetadata> chromosomes)
            {
                foreach (GenomeMetadata.SequenceMetadata referenceChromosome in chromosomes)
                {
                    if (chrom == referenceChromosome.Name && referenceChromosome.IsDecoyOrOther())
                    {
                        return true;
                    }
                }
                return false;
            }

            public bool IsDecoyOrOther()
            {
                switch (this.Type)
                {
                    case SequenceType.Decoy:
                    case SequenceType.Other:
                        return true;
                    default:
                        return false;
                }
            }

            /// <summary>
            /// Checks if the chromosome is mitochondrial
            /// Assumes species is Homo_sapiens 
            /// At some point this metadata should be part of the GenomeSize.xml
            /// </summary>
            /// <param name="chr">Name of chromosome to check.</param>
            /// <returns></returns>
            public static bool IsMito(string chr)
            {
                string tempChr = chr.ToLowerInvariant();
                if (tempChr == "chrm" || tempChr == "mt") return true;
                return false;
            }

            public bool IsMito()
            {
                switch (this.Type)
                {
                    case SequenceType.Mitochondria:
                        return true;
                    case SequenceType.Contig:
                    case SequenceType.Decoy:
                    case SequenceType.Autosome:
                    case SequenceType.Allosome:
                    case SequenceType.Other:
                        return false;
                    case SequenceType.Unknown:
                    default:
                        return IsMito(Name);
                }
            }

            /// <summary>
            ///     Used when sorting the sequences according to length (descending order)
            /// </summary>
            public int CompareTo(SequenceMetadata chromosome)
            {
                return chromosome.Length.CompareTo(Length);
            }

            /// <summary>
            ///     Checks the index file for the sequence offset
            /// </summary>
            public long GetFileOffset()
            {
                string faiPath = string.Format("{0}.fai", FastaPath);

                if (!File.Exists(faiPath))
                {
                    throw new ApplicationException(string.Format("Cannot open the FASTA index file ({0}) for reading.", faiPath));
                }

                long referenceOffset = 0;

                using (StreamReader faiReader = new StreamReader(faiPath))
                {
                    bool foundReference = false;

                    while (true)
                    {
                        // get the next line
                        string line = faiReader.ReadLine();
                        if (line == null) break;

                        // split the columns
                        string[] faiColumns = line.Split('\t');

                        if (faiColumns.Length != 5)
                        {
                            throw new ApplicationException(string.Format("Expected 5 columns in the FASTA index file ({0}), but found {1}.",
                                faiPath, faiColumns.Length));
                        }

                        // check the reference name
                        if (faiColumns[0] == Name)
                        {
                            referenceOffset = long.Parse(faiColumns[2]);
                            foundReference = true;
                            break;
                        }
                    }

                    // sanity check
                    if (!foundReference)
                    {
                        throw new ApplicationException(
                            string.Format("Unable to find the current sequence ({0}) in the index file ({1})", Name,
                                          faiPath));
                    }
                }

                return referenceOffset;
            }

            /// <summary>
            ///     Retrieves the bases associated with this sequence
            /// </summary>
            public string GetBases(bool toUpper = false)
            {
                long referenceOffset = GetFileOffset();
                long numRemainingBases = Length;

                StringBuilder builder = new StringBuilder();
                using (FileStream stream = new FileStream(FastaPath, FileMode.Open, FileAccess.Read, FileShare.Read, BufferSize))
                using (StreamReader reader = new StreamReader(stream))
                {
                    stream.Position = referenceOffset;
                    string line = string.Empty;

                    while (numRemainingBases > 0)
                    {
                        line = reader.ReadLine();
                        if (line == null)
                        {
                            throw new ApplicationException(string.Format(
                                    "Encountered the end of file before being able to retrieve the entire FASTA sequence. Remaining bases: {0}",
                                    numRemainingBases));
                        }
                        if (toUpper)
                        {
                            builder.Append(line.ToUpperInvariant());
                        }
                        else
                        {
                            builder.Append(line);
                        }
                        numRemainingBases -= line.Length;
                    }
                }

                return builder.ToString();
            }

            /// <summary>
            ///     Writes the current sequence to a FASTA file
            /// </summary>
            public void WriteFastaFile(string outputFastaPath)
            {
                using (FileStream readerFS = new FileStream(FastaPath, FileMode.Open, FileAccess.Read,
                    FileShare.Read, BufferSize))
                using (StreamReader reader = new StreamReader(readerFS))
                using (StreamWriter writer = new StreamWriter(new FileStream(outputFastaPath,
                    FileMode.Create, FileAccess.Write, FileShare.Write, BufferSize)))
                {
                    // initialize
                    string line = string.Empty;
                    writer.NewLine = "\n";
                    long numRemainingBases = Length;

                    // jump to the reference offset
                    readerFS.Position = GetFileOffset();

                    // write the header
                    writer.WriteLine(">{0}", Name);

                    // write the FASTA bases
                    while (numRemainingBases > 0)
                    {
                        line = reader.ReadLine();
                        if (line == null)
                        {
                            throw new ApplicationException(string.Format(
                                "Encountered the end of file before being able to retrieve the entire FASTA sequence. Remaining bases: {0}",
                                numRemainingBases));
                        }
                        writer.WriteLine(line);
                        numRemainingBases -= line.Length;
                    }
                }
            }
        }
    }

    /// <summary>
    ///     Stores additional information used by BOLT but not by ELAND
    /// </summary>
    public class FastaHeaderMetadata
    {
        public bool IsCircular = false;

        public FastaHeaderMetadata()
        {
            Ploidy = 2;
        }

        public int Ploidy { get; set; }
        public string Build { get; set; }
        public string Checksum { get; set; }
        public string Species { get; set; }
    }

    /// <summary>
    ///     Used to store the indices in the FASTA index (fai) file.
    ///     N.B. This is what samtools and BOLT uses as an index
    /// </summary>
    public class RefIndexEntry
    {
        public string Checksum;
        public string FilePath;
        public long KnownBasesLength;
        public long Length;
        public int LineLength;
        public string Name;
        public int NonWSLineLength;
        public long Offset;
    }

    /// <summary>
    ///     Parses a FASTA reference file and creates index (.fai), dictionary (.dict)
    /// </summary>
    public class ReferenceIndexer
    {
        #region Members
        private const int BufferSize = 10485760;
        private const string CircularTag = "CIRCULAR";
        private const byte NonCanonicalBase = 4;
        private readonly byte[] _convertBaseToNumber;
        private readonly byte[] _convertBaseToUpper;
        private readonly Regex _referenceBuildRegex = new Regex(@"BUILD\(([^\)]+)\)", RegexOptions.Compiled);
        private readonly Regex _referenceNameRegex = new Regex(@"^(\S+)", RegexOptions.Compiled);
        private readonly Regex _referencePloidyRegex = new Regex(@"PLOIDY\(([^\)]+)\)", RegexOptions.Compiled);
        private readonly Regex _referenceSpeciesRegex = new Regex(@"SPECIES\(([^\)]+)\)", RegexOptions.Compiled);
        private List<byte> _basesList = new List<byte>();
        private int _currentKnownBaseLineLen;
        private int _currentLineLength;
        private long _currentPosition;
        private string _fastaPath;
        private bool _inValidRegion;
        private bool _isFirstFastaEntry;
        private long _knownBaseRefLen;
        private HashSet<int> _lineLengths;
        private List<FastaHeaderMetadata> _metadataList;
        private long _numValidBases;
        private FileStream _reader;
        private List<RefIndexEntry> _referenceIndexes;
        private long _referenceLength;
        private int _referenceLineLength;
        private int _referenceLineLengthNoWhitespace;
        private string _referenceName = string.Empty;
        private long _referenceOffset;
        private uint _regionBegin;
        private uint _regionEnd;
        // global variable to keep the known bases (non-N) in ref seq in each group (that is reset with > in fasta file)
        #endregion

        public ReferenceIndexer()
        {
            // create our base to number lookup table
            _convertBaseToNumber = new byte[256];
            for (int i = 0; i < 256; i++) _convertBaseToNumber[i] = NonCanonicalBase;
            _convertBaseToNumber['a'] = 0;
            _convertBaseToNumber['A'] = 0;
            _convertBaseToNumber['c'] = 1;
            _convertBaseToNumber['C'] = 1;
            _convertBaseToNumber['g'] = 2;
            _convertBaseToNumber['G'] = 2;
            _convertBaseToNumber['t'] = 3;
            _convertBaseToNumber['T'] = 3;

            _convertBaseToUpper = new byte[256];
            for (int i = 0; i < 256; i++) _convertBaseToUpper[i] = (byte)i;
            _convertBaseToUpper['a'] = (byte)'A';
            _convertBaseToUpper['c'] = (byte)'C';
            _convertBaseToUpper['g'] = (byte)'G';
            _convertBaseToUpper['t'] = (byte)'T';
            _convertBaseToUpper['n'] = (byte)'N';
        }

        private void Initialize()
        {
            _isFirstFastaEntry = true;
            _inValidRegion = false;
            _regionBegin = 0;
            _regionEnd = 0;
            _currentPosition = 0;
            _numValidBases = 0;
            _referenceName = string.Empty;
            _referenceLength = 0;
            _referenceOffset = 0;
            _referenceLineLength = 0;
            _referenceLineLengthNoWhitespace = 0;
            _currentLineLength = 0;
            _basesList = new List<byte>();
            _lineLengths = new HashSet<int>();
            _knownBaseRefLen = 0;
            _currentKnownBaseLineLen = 0;
        }

        /// <summary>
        ///     Adds the MD5 checksum to our FASTA header metadata list
        /// </summary>
        private static void AddMD5Checksum(ref List<byte> basesList, List<FastaHeaderMetadata> metadataList)
        {
            // compute the checksum
            MD5CryptoServiceProvider md5Provider = new MD5CryptoServiceProvider();
            byte[] data = basesList.ToArray();
            basesList.Clear();
            data = md5Provider.ComputeHash(data);

            // add the checksum to the last metadata entry
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < data.Length; i++) sb.Append(data[i].ToString("x2").ToLower());
            metadataList[metadataList.Count - 1].Checksum = sb.ToString();
        }

        private static string ReadLine(FileStream fs)
        {
            StringBuilder sb = new StringBuilder();

            // keep reading until we hit a carriage return
            while (true)
            {
                int i = fs.ReadByte();
                if (i == -1)
                {
                    throw new ApplicationException(
                        string.Format("EOF encountered before the end of the FASTA header in {0}", fs.Name));
                }
                if (IsEOL(i)) break;
                sb.Append((char)i);
            }

            return sb.ToString();
        }

        private static void UpdateValidRegions(ref bool inValidRegion, ref long numValidBases,
            long currentPosition, uint regionBegin)
        {
            if (inValidRegion)
            {
                uint regionEnd = (uint)(currentPosition - 1);
                numValidBases += regionEnd - regionBegin + 1;
                inValidRegion = false;
            }
        }

        private static void AddReferenceIndexEntry(ref List<RefIndexEntry> indices, string name, long length,
                                                   long offset, int lineLength)
        {
            RefIndexEntry ri = new RefIndexEntry
            {
                Length = length,
                LineLength = lineLength + 1,
                Name = name,
                NonWSLineLength = lineLength,
                Offset = offset
            };
            indices.Add(ri);
        }

        private void HandleRecordStart()
        {
            // handle the FASTA header
            if (!_isFirstFastaEntry)
            {
                AddMD5Checksum(ref _basesList, _metadataList);
                UpdateValidRegions(ref _inValidRegion, ref _numValidBases, _currentPosition, _regionBegin);

                if (_lineLengths.Count > 2)
                {
                    throw new ApplicationException(string.Format("Found inconsistent line lengths in reference {0} in file {1}", _referenceName,
                        _fastaPath));
                }

                // add the reference index information
                RefIndexEntry indexEntry = new RefIndexEntry
                {
                    Length = _referenceLength,
                    LineLength = _referenceLineLength,
                    Name = _referenceName,
                    NonWSLineLength = _referenceLineLengthNoWhitespace,
                    Offset = _referenceOffset,
                    Checksum = _metadataList[_metadataList.Count - 1].Checksum,
                    FilePath = _fastaPath,
                    KnownBasesLength = _knownBaseRefLen
                };
                _referenceIndexes.Add(indexEntry);
            }

            // look for empty headers
            string header = ReadLine(_reader);
            if (string.IsNullOrEmpty(header))
            {
                throw new ApplicationException(string.Format("Empty FASTA header found in {0}", _fastaPath));
            }

            // extract the reference name
            Match headerMatch = _referenceNameRegex.Match(header);
            if (!headerMatch.Success)
            {
                throw new ApplicationException(
                    string.Format("Unable to extract the reference name in {0} from line: {1}", _fastaPath, header));
            }

            _referenceName = headerMatch.Groups[1].Value;
            _referenceOffset = -1;
            _referenceLength = 0;
            _knownBaseRefLen = 0; // reset
            _lineLengths.Clear();

            // extract our metadata
            FastaHeaderMetadata fhm = new FastaHeaderMetadata();

            Match referenceBuildMatch = _referenceBuildRegex.Match(header);
            Match referenceSpeciesMatch = _referenceSpeciesRegex.Match(header);
            Match referencePloidyMatch = _referencePloidyRegex.Match(header);

            if (referenceBuildMatch.Success) fhm.Build = referenceBuildMatch.Groups[1].Value;
            if (referenceSpeciesMatch.Success) fhm.Species = referenceSpeciesMatch.Groups[1].Value;
            if (referencePloidyMatch.Success) fhm.Ploidy = int.Parse(referencePloidyMatch.Groups[1].Value);

            if (header.IndexOf(CircularTag) != -1) fhm.IsCircular = true;
            _metadataList.Add(fhm);

            if (_isFirstFastaEntry) _isFirstFastaEntry = false;
        }

        public void IndexFASTAFile(string inputFastaPath, List<FastaHeaderMetadata> metadata,
            List<RefIndexEntry> refIndexes)
        {
            // check that our FASTA file exists
            _fastaPath = Path.GetFullPath(inputFastaPath);
            if (!File.Exists(_fastaPath))
            {
                throw new ApplicationException(string.Format("Cannot open the FASTA file ({0}) for reading.", _fastaPath));
            }

            // derive our output filenames
            string outputDirectory = Path.GetDirectoryName(_fastaPath);
            _metadataList = metadata;
            _referenceIndexes = refIndexes;
            Initialize();

            try
            {
                using (_reader = new FileStream(_fastaPath, FileMode.Open, FileAccess.Read, FileShare.Read, BufferSize))
                {
                    int b = 0;
                    while (true)
                    {
                        // grab the next byte
                        int previousChar = b;
                        b = _reader.ReadByte();
                        if (b == -1) break; // End of file

                        if (IsEOL(b))
                        {
                            if (_currentLineLength > 0)
                            {
                                if (_referenceLength == 0)
                                {
                                    _referenceLineLength = _currentLineLength + 1;
                                    _referenceLineLengthNoWhitespace = _currentLineLength;
                                }
                                _lineLengths.Add(_currentLineLength);
                                _referenceLength += _currentLineLength;
                                _knownBaseRefLen += _currentKnownBaseLineLen;
                                _currentLineLength = 0;
                                _currentKnownBaseLineLen = 0; // reset;
                            }
                            else if (IsEOL(previousChar) && _referenceLength == _referenceLineLengthNoWhitespace)
                            {
                                // Second EOL character, on the first line of the reference:
                                _referenceLineLength++;
                            }
                            continue;
                        }
                        if (b == '>')
                        {
                            HandleRecordStart();
                        }
                        else
                        {
                            if (_referenceOffset == -1) _referenceOffset = _reader.Position - 1;
                            // handle the bases
                            _currentLineLength++;

                            // handle the unknown base (N)
                            // this is not elegant like Reg expresion, but it's faster
                            switch (b)
                            {
                                case 'A':
                                case 'a':
                                case 'C':
                                case 'c':
                                case 'G':
                                case 'g':
                                case 'T':
                                case 't':
                                    _currentKnownBaseLineLen++;
                                    break;
                            }

                            uint num = _convertBaseToNumber[b];
                            _basesList.Add(_convertBaseToUpper[b]);

                            if (num == NonCanonicalBase)
                            {
                                num = 0; // A

                                if (_inValidRegion)
                                {
                                    _regionEnd = (uint)(_currentPosition - 1);
                                    _numValidBases += _regionEnd - _regionBegin + 1;
                                    _inValidRegion = false;
                                }
                            }
                            else
                            {
                                if (!_inValidRegion)
                                {
                                    _regionBegin = (uint)_currentPosition;
                                    _inValidRegion = true;
                                }
                            }

                            _currentPosition++;

                        }
                    } // Loop over bytes b


                    // write the last valid region
                    if (_currentLineLength > 0)
                    {
                        if (_referenceLength == 0)
                        {
                            // This case triggers for a single-line reference (all bases on same line, so no trailing whitespace):
                            _referenceLineLength = _currentLineLength;
                            _referenceLineLengthNoWhitespace = _currentLineLength;
                        }
                        _lineLengths.Add(_currentLineLength);
                        _referenceLength += _currentLineLength;
                        _knownBaseRefLen += _currentKnownBaseLineLen;
                    }

                    if (_lineLengths.Count > 2)
                    {
                        throw new ApplicationException(string.Format(
                            "Found inconsistent line lengths in reference {0} in file {1}", _referenceName, _fastaPath));
                    }

                    AddMD5Checksum(ref _basesList, _metadataList);
                    UpdateValidRegions(ref _inValidRegion, ref _numValidBases, _currentPosition, _regionBegin);

                    // add the reference index information
                    RefIndexEntry rie = new RefIndexEntry
                    {
                        Length = _referenceLength,
                        LineLength = _referenceLineLength,
                        Name = _referenceName,
                        NonWSLineLength = _referenceLineLengthNoWhitespace,
                        Offset = _referenceOffset,
                        Checksum = _metadataList[_metadataList.Count - 1].Checksum,
                        FilePath = _fastaPath,
                        KnownBasesLength = _knownBaseRefLen
                    };
                    _referenceIndexes.Add(rie);
                }

                // create the FASTA index (fai) and GATK dict files
                WriteDictFile(outputDirectory);
                WriteIndexFile();
            }
            catch (Exception e)
            {
                throw new ApplicationException(
                    string.Format("An exception occurred when packing the reference sequences: {0}", e.Message));
            }
        }

        private void WriteIndexFile()
        {
            string faiFilename = _fastaPath + ".fai";
            using (StreamWriter faiWriter = new StreamWriter(new FileStream(faiFilename, FileMode.Create)))
            {
                faiWriter.NewLine = "\n";
                foreach (RefIndexEntry ri in _referenceIndexes)
                {
                    faiWriter.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}",
                        ri.Name, ri.Length, ri.Offset, ri.NonWSLineLength, ri.LineLength);
                }
            }
        }

        private void WriteDictFile(string outputDirectory)
        {
            string dictFilename = Path.Combine(outputDirectory, Path.GetFileNameWithoutExtension(_fastaPath)) + ".dict";
            using (StreamWriter dictWriter = new StreamWriter(new FileStream(dictFilename, FileMode.Create)))
            {
                dictWriter.NewLine = "\n";
                dictWriter.WriteLine("@HD\tVN:1.0\tSO:unsorted");
                foreach (RefIndexEntry ri in _referenceIndexes)
                {
                    dictWriter.WriteLine("@SQ\tSN:{0}\tLN:{1}\tUR:file:{2}\tM5:{3}",
                        ri.Name, ri.Length, ri.FilePath, ri.Checksum);
                }
            }
        }

        /// <summary>
        ///     Create .fai and .dict files for the specified FASTA reference file.  
        /// </summary>
        public static void CreateFASTAIndexFiles(string fastaPath, List<FastaHeaderMetadata> metadataList,
            List<RefIndexEntry> referenceIndexes)
        {
            ReferenceIndexer indexer = new ReferenceIndexer();
            indexer.IndexFASTAFile(fastaPath, metadataList, referenceIndexes);
        }

        /// <summary>
        ///     Detects if a character represents the end of the line
        /// </summary>
        /// <returns>true if the character is a carriage return or line feed</returns>
        public static bool IsEOL(int b)
        {
            if ((b == '\n') || (b == '\r')) return true;
            return false;
        }
    }
}