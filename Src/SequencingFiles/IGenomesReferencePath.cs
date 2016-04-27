using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace SequencingFiles
{
    /// <summary>
    /// This class provides access to the Species, Provider and Build from a path to the reference genome
    /// Assume reference path is in iGenomes format: 
    ///   "...\{species}\{provider}\{build}\Sequence\WholeGenomeFasta" or
    ///   "...\{species}\{provider}\{build}\Sequence\WholeGenomeFasta\genome.fa"
    /// This class does not impose a restriction on the species name, but
    /// typically it follows binomial nomenclature with an underscore between the two names (e.g. "Homo_species") 
    /// 
    /// Developer Note:
    ///		If you have access to the GenomeMetadata it is recommended to use that class to get at the Species and Build metadata
    ///		Currently, the GenomeMetadata class uses this class internally to determine that information
    ///		
    /// Example Use Case from the T/N workflow:
    ///		We have two bam files with headers that contain iGenomes style path to the reference genome.
    ///		These bam files were potentially generated on different systems with distinct paths (e.g. BaseSpace vs local cluster) 
    ///		We want to validate that both reference genomes refer to the same build before we continue with our analysis
    /// 			
    ///			string normalFasta = GetFastaFromBamHeader(normalSample.AlignmentPath);
    ///			string tumorFasta = GetFastaFromBamHeader(tumorSample.AlignmentPath);
    ///			IGenomesReferencePath normalReference = IGenomesReferencePath.GetReferenceFromFastaPath(normalFasta);
    ///			IGenomesReferencePath tumorReference = IGenomesReferencePath.GetReferenceFromFastaPath(tumorFasta);
    ///			if (!normalReference.Equals(tumorReference))
    ///				throw new ApplicationException(string.Format("Incompatible reference genomes in the tumor and normal bam headers. Analysis cannot proceed. Paths were {0} and {1}", normalFasta, tumorFasta));
    /// 
    /// </summary>
    public class IGenomesReferencePath
    {
        private IGenomesReferencePath(string referencePath)
        {
            if (string.IsNullOrEmpty(referencePath))
                throw new ArgumentException("Reference path cannot be null or empty");

            if (referencePath.StartsWith("file:"))
            {
                try
                {
                    referencePath = Path.GetFullPath(new Uri(referencePath, UriKind.RelativeOrAbsolute).LocalPath).TrimEnd(new[] { '/', '\\' });
                }
                catch (Exception)
                {
                    referencePath = Path.GetFullPath(referencePath.TrimEnd(new[] { '/', '\\' }));
                }
            }
            else
            {
                referencePath = Path.GetFullPath(referencePath.TrimEnd(new[] { '/', '\\' }));
            }

            int offset = 5;
            if (referencePath.EndsWith(".fa"))
                offset++;

            string[] Bits = referencePath.Split(Path.DirectorySeparatorChar);
            bool hasSequence = false;
            bool hasWholeGenomeFasta = false;
            if (Bits.Length >= offset)
            {
                hasSequence = Bits[Bits.Length - offset + 3].ToLowerInvariant() == "sequence";
                hasWholeGenomeFasta = Bits[Bits.Length - offset + 4].ToLowerInvariant() == "wholegenomefasta";
                // Note: BWA header may point to Sequence/BWAIndex/foo.fa instead of Sequence/WholeGenomeFasta/foo.fa.  That's okay:
                if (!hasWholeGenomeFasta && Bits[Bits.Length - offset + 4].ToLowerInvariant() == "bwaindex")
                    hasWholeGenomeFasta = true;
            }
            if (hasSequence && hasWholeGenomeFasta)
            {
                Species = Bits[Bits.Length - offset];
                Provider = Bits[Bits.Length - offset + 1];
                Build = Bits[Bits.Length - offset + 2];
            }
            else
            {
            throw new ArgumentException("Reference path not in iGenomes format ...\\{species}\\{provider}\\{version}\\Sequence\\WholeGenomeFasta or ...\\{species}\\{provider}\\{version}\\Sequence\\WholeGenomeFasta\\genome.fa." + string.Format("Input was: {0}", referencePath));
            }

        }

        public string Species { get; private set; }
        public string Provider { get; private set; }
        public string Build { get; private set; }

        private static IGenomesReferencePath SafeGetReference(string fastaPath)
        {
            IGenomesReferencePath reference = null;
            try
            {
                reference = new IGenomesReferencePath(fastaPath);
            }
            catch { }
            return reference;
        }

        public static IGenomesReferencePath GetReferenceFromFastaPath(string fastaPath)
        {
            return SafeGetReference(fastaPath);
        }

        public static IGenomesReferencePath GetReference(SequencingFiles.GenomeMetadata genome)
        {
            return GetReference(genome.Sequences[0]);
        }

        public static IGenomesReferencePath GetReference(SequencingFiles.GenomeMetadata.SequenceMetadata sequence)
        {
            return SafeGetReference(sequence.FastaPath);
        }

        public static IGenomesReferencePath GetReferenceFromVcfHeader(string vcfPath)
        {
            string referencePath = null;
            using (GzipReader Reader = new GzipReader(vcfPath))
            {
                while (true)
                {
                    string FileLine = Reader.ReadLine();
                    if (FileLine == null || !FileLine.StartsWith("#"))
                        break;
                    if (FileLine.StartsWith("##reference="))
                        referencePath = FileLine.Substring(12);
                }
            }
            return SafeGetReference(referencePath);
        }

        public static IGenomesReferencePath GetReferenceFromBamHeader(string bamPath)
        {
            using (BamReader reader = new BamReader(bamPath))
            {
                return SafeGetReference(reader.GetReferenceURI());
            }
        }

        public bool Equals(IGenomesReferencePath p)
        {
            return p.Species == Species && p.Provider == Provider && p.Build == Build;
        }

        public string GetWholeGenomeFasta(string genomesDirectory)
        {
            if (Species == null || Provider == null || Build == null)
                return null;
            string genomePath = Path.Combine(genomesDirectory, Species, Provider, Build, "Sequence", "WholeGenomeFasta", "genome.fa");
            return genomePath;
        }

        public override string ToString()
        {
            return string.Format("{0} ({1} {2})", Species.Replace("_", " "), Provider, Build);
        }

        public static bool SameReference(string fastaA, string fastaB)
        {
            IGenomesReferencePath fastaAReference = null;
            IGenomesReferencePath fastaBReference = null;
            try
            {
                fastaAReference = new IGenomesReferencePath(fastaA);
            }
            catch (ArgumentException)
            { }
            try
            {
                fastaBReference = new IGenomesReferencePath(fastaB);
            }
            catch (ArgumentException)
            { }

            if (fastaAReference == null || fastaBReference == null)
                return fastaA == fastaB;
            return fastaAReference.Equals(fastaBReference);
        }
    }
}
