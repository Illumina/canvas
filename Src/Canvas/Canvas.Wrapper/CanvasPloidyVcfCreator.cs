using Illumina.Common.FileSystem;
using Illumina.SecondaryAnalysis.VariantCalling;
using Isas.Framework.DataTypes;
using Isas.SequencingFiles;
using System.Collections.Generic;
using Isas.Ploidy;
using System;
using System.Linq;

namespace Canvas.Wrapper
{
    public class CanvasPloidyVcfCreator
    {
        public const string PloidyVcfName = "ploidy.vcf.gz";
        public CanvasPloidyVcfCreator(PloidyCorrector ploidyFixer)
        {
            _ploidyFixer = ploidyFixer;
        }
        private readonly PloidyCorrector _ploidyFixer;

        /// <summary>
        /// Write out the ploidy vcf file if ploidy information is available from the vcf header
        /// </summary>
        public Vcf CreatePloidyVcf(SampleSet<SexPloidyInfo> ploidyInfos, GenomeMetadata genomeMetadata, IDirectoryLocation sampleSandbox)
        {
            var ploidyVcf = new Vcf(sampleSandbox.GetFileLocation(PloidyVcfName));
            _ploidyFixer.WritePloidyVcfFile(ploidyVcf, ploidyInfos, genomeMetadata);
            return ploidyVcf;
        }

        public Vcf CreatePloidyVcf(string sampleId, SexPloidyInfo sexPloidyInfo, GenomeMetadata genomeMetadata, IDirectoryLocation sampleSandbox)
        {
            var sampleInfo = new SampleInfo(sampleId, "SampleName");
            var sampleSet = new SampleSet<SexPloidyInfo>(new Dictionary<SampleInfo, SexPloidyInfo> { { sampleInfo, sexPloidyInfo } });
            return CreatePloidyVcf(sampleSet, genomeMetadata, sampleSandbox);
        }

        public Vcf CreatePloidyVcfFromSexChromosomeKaryotype(string sampleId, string SexChromosomeKaryotype, GenomeMetadata genomeMetadata, IDirectoryLocation sampleSandbox)
        {

            SexPloidyInfo sexPloidyInfo = GetSexPloidyFromSexChromosomeKaryotype(SexChromosomeKaryotype);
            return CreatePloidyVcf(sampleId, sexPloidyInfo, genomeMetadata, sampleSandbox);
        }

        private static SexPloidyInfo GetSexPloidyFromSexChromosomeKaryotype(string sexChromosomeKaryotype)
        {
            var karyotypeLower = sexChromosomeKaryotype.ToLowerInvariant();
            var ploidyX = karyotypeLower.Count(letter => letter == 'x');
            var ploidyY = karyotypeLower.Count(letter => letter == 'y');
            return new SexPloidyInfo(ploidyX, ploidyY);
        }
    }
}