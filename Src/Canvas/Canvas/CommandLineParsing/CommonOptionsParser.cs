using Canvas.CommandLineParsing.CoreOptionTypes;
using Canvas.CommandLineParsing.OptionProcessing;

namespace Canvas.CommandLineParsing
{
    public class CommonOptionsParser : Option<CommonOptions>
    {
        private static readonly FileOption BAlleleSites = FileOption.CreateRequired("vcf containing SNV b-allele sites (only sites with PASS in the filter column will be used)", "b-allele-vcf");
        private static readonly FlagOption IsDbSnpVcf = new FlagOption("exclude SNV b-allele sites from the vcf that do not have sufficient read evidence for a heterozygous genotype in the sample. " +
                                                                       "This option should be used when the b-allele vcf is not matched to the sample (e.g. from dbSNP)", "exclude-non-het-b-allele-sites");
        private static readonly FileOption PloidyBed = FileOption.Create(".bed file containing regions of known ploidy in the sample. Copy number calls matching the known ploidy in these regions will be considered non-variant", "ploidy-bed");
        private static readonly FileOption KmerFasta = FileOption.CreateRequired("Canvas-ready reference fasta file", "r", "reference");
        private static readonly DirectoryOption Output = DirectoryOption.CreateRequired("output directory", "o", "output");
        private static readonly DirectoryOption WholeGenomeFasta = DirectoryOption.CreateRequired("folder that contains both genome.fa and GenomeSize.xml", "g", "genome-folder");
        private static readonly FileOption FilterFile = FileOption.CreateRequired(".bed file of regions to skip", "f", "filter-bed");
        private static readonly StringOption SampleName = StringOption.CreateRequired("sample name", "n", "sample-name");
        private static readonly DictionaryOption CustomParameters = DictionaryOption.Create("use custom parameter for command-line tool. VALUE must contain the tool name followed by a comma and then the custom parameters.", "custom-parameters");
        private static readonly StringOption StartCheckpoint = StringOption.Create("continue analysis starting at the specified checkpoint. VALUE can be the checkpoint name or number", "c");
        private static readonly StringOption StopCheckpoint = StringOption.Create("stop analysis after the specified checkpoint is complete. VALUE can be the checkpoint name or number", "s");

        public override OptionCollection<CommonOptions> GetOptions()
        {
            return new OptionCollection<CommonOptions>
            {
                BAlleleSites, IsDbSnpVcf, PloidyBed, Output, KmerFasta, WholeGenomeFasta, FilterFile, SampleName, CustomParameters, //StartCheckpoint, StopCheckpoint
            };
        }

        public override ParsingResult<CommonOptions> Parse(SuccessfulResultCollection parseInput)
        {
            var bAlleleSites = parseInput.Get(BAlleleSites);
            bool isDbSnpVcf = parseInput.Get(IsDbSnpVcf);
            var ploidyBed = parseInput.Get(PloidyBed);
            var output = parseInput.Get(Output);
            var wholeGenomeFasta = parseInput.Get(WholeGenomeFasta);
            var kmerFasta = parseInput.Get(KmerFasta);
            var filterBed = parseInput.Get(FilterFile);
            var sampleName = parseInput.Get(SampleName);
            var customParameters = parseInput.Get(CustomParameters);
            string startCheckpoint = null;//parseInput.Get(StartCheckpoint);
            string stopCheckpoint = null;//parseInput.Get(StopCheckpoint);
            return ParsingResult<CommonOptions>.SuccesfulResult(
                new CommonOptions(
                    bAlleleSites,
                    isDbSnpVcf,
                    ploidyBed,
                    output,
                    wholeGenomeFasta,
                    kmerFasta,
                    filterBed,
                    sampleName,
                    customParameters,
                    startCheckpoint,
                    stopCheckpoint));
        }
    }
}