using System;

namespace Isas.Shared
{
    public interface IGtfAnnotationResource
    {
        IFileLocation GTFFile { get; }
    }

    public interface IRnaReferenceGenome : IReferenceGenome, IGtfAnnotationResource
    {
        IDirectoryLocation GenesFolder { get; }
    }

    public class RnaReferenceGenome : ReferenceGenome, IRnaReferenceGenome
    {
        private bool _validate;

        public IDirectoryLocation GenesFolder { get; }

        public RnaReferenceGenome(IDirectoryLocation genomeBuildDirectory,
            IDirectoryLocation genesDir, bool validate = true)
            : base(genomeBuildDirectory)
        {
            GenesFolder = genesDir;

            _validate = validate;
        }

        public IFileLocation GTFFile
        {
            get
            {
                IFileLocation f = GenesFolder.GetFileLocation("genes.gtf");
                if (_validate && !f.Exists)
                    throw new ApplicationException("Missing annotation file: " + f.FullName);
                return f;
            }
        }

        // These methods to be moved to STAR wrappers
        public IFileLocation AugmentedFastaFile
        {
            get
            {
                IFileLocation f = SequenceDirectory.GetFileLocation("STARSequence", "genome.augmented.fa");
                if (_validate && !f.Exists)
                    throw new ApplicationException("Missing genome sequence file: " + f.FullName);
                return f;
            }
        }
        public IFileLocation ChromList
        {
            get
            {
                IFileLocation f = SequenceDirectory.GetFileLocation("STARSequence", "chromList.txt");
                return (_validate && !f.Exists) ? null : f;

            }
        }
    }
}
