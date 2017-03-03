using CanvasCommon.CommandLineParsing.CoreOptionTypes;
using CanvasCommon.CommandLineParsing.OptionProcessing;

namespace Canvas.CommandLineParsing
{
    internal class SingleSampleCommonOptionsParser : Option<SingleSampleCommonOptions>
    {
        public const string SampleBAlleleVcfOptionName = "sample-b-allele-vcf";
        public const string PopulationBAlleleVcfOptionName = "population-b-allele-vcf";
        public const string PloidyBedOptionName = "ploidy-bed";
        private static readonly FileOption SampleBAlleleSites = FileOption.CreateRequired("vcf containing SNV b-allele sites in the sample (only sites with PASS in the filter column will be used)", SampleBAlleleVcfOptionName);
        public static readonly FileOption PopulationBAlleleSites = FileOption.CreateRequired("vcf containing SNV b-allele sites in the population (only sites with PASS in the filter column will be used)", PopulationBAlleleVcfOptionName);
        private static readonly ExclusiveFileOption BAlleleSites = ExclusiveFileOption.CreateRequired(SampleBAlleleSites, PopulationBAlleleSites);
        private static readonly FileOption PloidyBed = FileOption.Create(".bed file containing regions of known ploidy in the sample. Copy number calls matching the known ploidy in these regions will be considered non-variant", PloidyBedOptionName);
        private static readonly StringOption SampleName = StringOption.CreateRequired("sample name", "n", "sample-name");

        public override OptionCollection<SingleSampleCommonOptions> GetOptions()
        {
            return new OptionCollection<SingleSampleCommonOptions>
            {
                BAlleleSites, PloidyBed, SampleName
            };
        }

        public override ParsingResult<SingleSampleCommonOptions> Parse(SuccessfulResultCollection parseInput)
        {
            var bAlleleSites = parseInput.Get(BAlleleSites);
            bool isDbSnpVcf = bAlleleSites.MatchedOption.Equals(PopulationBAlleleSites);
            var ploidyBed = parseInput.Get(PloidyBed);
            var sampleName = parseInput.Get(SampleName);
            return ParsingResult<SingleSampleCommonOptions>.SuccessfulResult(
                new SingleSampleCommonOptions(
                    bAlleleSites.Result,
                    isDbSnpVcf,
                    ploidyBed,
                    sampleName));
        }
    }
}