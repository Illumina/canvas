using System;
using SequencingFiles;
using System.IO;

namespace Isas.Shared.Utilities
{
    /// <summary>
    /// Provides methods for navigating the IGenomes directory tree
    /// </summary>
    public interface IReferenceGenome
    {
        string Species { get; }
        string Provider { get; }
        string Build { get; }

        IFileLocation WholeGenomeFasta { get; }
        IFileLocation WholeGenomeFastaIndex { get; }
        GenomeMetadata GenomeMetadata { get; }

        IDirectoryLocation AnnotationDirectory { get; }
        IDirectoryLocation SequenceDirectory { get; }
        IDirectoryLocation WholeGenomeFastaDirectory { get; }
        IFileLocation GenomeSizeXml { get; }
    }

    public interface IReferenceGenomeFactory
    {
        IReferenceGenome GetReferenceGenome(IDirectoryLocation wholeGenomeFasta);
    }

    public class ReferenceGenomeFactory : IReferenceGenomeFactory
    {
        public IReferenceGenome GetReferenceGenome(IDirectoryLocation wholeGenomeFasta)
        {
            return new ReferenceGenome(wholeGenomeFasta.Parent.Parent);
        }
    }

    public class ReferenceGenome : IReferenceGenome
    {
        public string Species => BuildDirectory.Parent.Parent.Name;
        public string Provider => BuildDirectory.Parent.Name;
        public string Build => BuildDirectory.Name;
        private IDirectoryLocation BuildDirectory { get; }
        public IDirectoryLocation AnnotationDirectory => BuildDirectory.GetDirectoryLocation("Annotation");
        public IDirectoryLocation SequenceDirectory => BuildDirectory.GetDirectoryLocation("Sequence");
        public IDirectoryLocation WholeGenomeFastaDirectory => SequenceDirectory.GetDirectoryLocation("WholeGenomeFasta");
        public IFileLocation WholeGenomeFasta => WholeGenomeFastaDirectory.GetFileLocation("genome.fa");
        public IFileLocation WholeGenomeFastaIndex => WholeGenomeFastaDirectory.GetFileLocation("genome.fai");
        public IFileLocation GenomeSizeXml => WholeGenomeFastaDirectory.GetFileLocation("GenomeSize.xml");

        public GenomeMetadata GenomeMetadata
        {
            get
            {
                var genomeMetadata = new GenomeMetadata();
                genomeMetadata.Deserialize(GenomeSizeXml.FullName);
                return genomeMetadata;
            }
        }

        public IFileLocation GetAnnotation(params string[] args)
        {
            return AnnotationDirectory.GetFileLocation(args);
        }

        public IFileLocation GetSequence(params string[] args)
        {
            return SequenceDirectory.GetFileLocation(args);
        }

        public ReferenceGenome(IDirectoryLocation buildDirectory)
        {
            BuildDirectory = buildDirectory;
        }

        public override string ToString()
        {
            return $"{Species}/{Provider}/{Build}";
        }
    }
}
