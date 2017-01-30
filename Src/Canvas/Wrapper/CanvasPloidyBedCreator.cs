namespace Illumina.SecondaryAnalysis.VariantCalling.StructuralVariants.Canvas
{
    public class CanvasPloidyBedCreator
    {
        public CanvasPloidyBedCreator(ILogger logger, IWorkManager workManager, PloidyCorrector ploidyFixer)
        {
            _logger = logger;
            _workManager = workManager;
            _ploidyFixer = ploidyFixer;
        }

        private readonly IWorkManager _workManager;
        private readonly ILogger _logger;
        private readonly PloidyCorrector _ploidyFixer;

        /// <summary>
        /// Write out the ploidy bed file if ploidy information is available from the vcf header
        /// </summary>
        public IFileLocation CreatePloidyBed(Vcf vcf, GenomeMetadata genomeMetadata, IDirectoryLocation sampleSandbox)
        {
            IFileLocation ploidyBed = sampleSandbox.GetFileLocation("ploidy.bed.gz");
            string fastaPath = genomeMetadata.Sequences.First().FastaPath;
            if (_ploidyFixer.GeneratePloidyBedFileFromVcf(
                genomeMetadata,
                fastaPath,
                vcf.VcfFile.FullName,
                ploidyBed.FullName, sampleSandbox.FullName, _logger, _workManager))
            {
                return ploidyBed;
            }
            _logger.Warn($"Sex chromosome ploidy not found in {vcf.VcfFile} header. No ploidy will be provided to Canvas.");
            return null;
        }

        /// <summary>
        ///  Write out the ploidy bed file if ploidy information is available from the vcf header
        /// Only create the normal XX or XY ploidy bed file so that Canvas can properly classify any abnormalities as variant. 
        /// If ploidy Y is > 1 produce the XY ploidy bed file, otherwise produce the XX ploidy bed file
        /// </summary>
        public IFileLocation CreateGermlinePloidyBed(Vcf vcf, GenomeMetadata genomeMetadata, IDirectoryLocation sampleSandbox)
        {
            string sexKaryotype = PloidyCorrector.GetSexChromosomeKaryotypeFromVcfHeader(vcf.VcfFile.FullName);
            if (sexKaryotype == null)
            {
                _logger.Warn($"Sex chromosome ploidy not found in {vcf.VcfFile} header. No ploidy will be provided to Canvas.");
                return null;
            }
            _logger.Info($"Found sex chromosome ploidy {PloidyCorrector.PrintPloidy(sexKaryotype)} in {vcf.VcfFile}");
            var ploidyInfo = new SamplePloidyInfo();
            IFileLocation ploidyBed = sampleSandbox.GetFileLocation("ploidy.bed.gz");
            if (sexKaryotype.ToLower().Contains("y"))
            {
                ploidyInfo.ProvidedPloidyX = 1;
                ploidyInfo.ProvidedPloidyY = 1;
                _logger.Info($"Creating male ploidy bed file at {ploidyBed}.");
            }
            else
            {
                ploidyInfo.ProvidedPloidyX = 2;
                ploidyInfo.ProvidedPloidyY = 0;
                _logger.Info($"Creating female ploidy bed file at {ploidyBed}.");
            }
            string headerLine = $"{PloidyCorrector.ReferenceSexChromosomeKaryotype}={PloidyCorrector.PrettyPrintPloidy(ploidyInfo.ProvidedPloidyX.Value, ploidyInfo.ProvidedPloidyY.Value)}";
            _ploidyFixer.WritePloidyBedFile(ploidyInfo, genomeMetadata, _ploidyFixer.GetParRegions(genomeMetadata),
                ploidyBed.FullName, headerLine, ploidy => true);
            return ploidyBed;
        }

        public bool GeneratePloidyBedFileFromSexChromosomeKaryotype(GenomeMetadata genomeMetadata, string fastaPath, string sexChromosomeKaryotype, string ploidyBedPath, string logFolder)
        {
            return _ploidyFixer.GeneratePloidyBedFileFromSexChromosomeKaryotype(genomeMetadata, fastaPath, sexChromosomeKaryotype, ploidyBedPath, logFolder);
        }
    }
}