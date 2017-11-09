using System;
using System.IO;
using System.Linq;
using Illumina.Common.FileSystem;
using Isas.SequencingFiles;

namespace Canvas.Wrapper
{
    public interface ICanvasAnnotationFileProvider
    {
        bool IsSupported(GenomeMetadata genome);
        bool CustomDbSnpVcf(GenomeMetadata genome);
        IFileLocation GetDbSnpVcf(GenomeMetadata genome);
        IFileLocation GetKmerFasta(GenomeMetadata genome);

        IFileLocation GetFilterBed(GenomeMetadata genome);

        IFileLocation GetCanvasAnnotationFile(GenomeMetadata genome, string fileName);
    }

    public class CanvasAnnotationFileProvider : ICanvasAnnotationFileProvider
    {
        private readonly IFileLocation _dbSnpVcf;
        private readonly IReferenceGenomeFactory _referenceGenomeFactory;

        public CanvasAnnotationFileProvider(IFileLocation dbSnpVcf, IReferenceGenomeFactory referenceGenomeFactory)
        {
            _dbSnpVcf = dbSnpVcf;
            _referenceGenomeFactory = referenceGenomeFactory;
        }

        public bool IsSupported(GenomeMetadata genome)
        {
            string species = genome.Contigs().FirstOrDefault()?.Species;
            return "Homo_sapiens".Equals(species, StringComparison.OrdinalIgnoreCase);
        }

        public bool CustomDbSnpVcf(GenomeMetadata genome)
        {
            return _dbSnpVcf != null;
        }

        public IFileLocation GetDbSnpVcf(GenomeMetadata genome)
        {
            if (_dbSnpVcf != null) return _dbSnpVcf;

            return GetCanvasAnnotationFile(genome, "dbsnp.vcf");
        }

        public IFileLocation GetKmerFasta(GenomeMetadata genome)
        {
            return GetCanvasAnnotationFile(genome, "kmerv2.fa");
        }

        public IFileLocation GetFilterBed(GenomeMetadata genome)
        {
            return GetCanvasAnnotationFile(genome, "filter13.bed");
        }

        public IFileLocation GetCanvasAnnotationFile(GenomeMetadata genome, string fileName)
        {
            return GetReferenceGenomeFromGenomeMetadata(genome).AnnotationDirectory.GetFileLocation("Canvas", fileName);
        }

        private IReferenceGenome GetReferenceGenomeFromGenomeMetadata(GenomeMetadata genomeMetadata)
        {
            return _referenceGenomeFactory.GetReferenceGenome(new DirectoryLocation(Path.GetDirectoryName(genomeMetadata.Contigs().First().FastaPath)));
        }
    }
}
