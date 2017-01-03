using CanvasCommon.CommandLineParsing.CoreOptionTypes;
using CanvasCommon.CommandLineParsing.OptionProcessing;

namespace Canvas.CommandLineParsing
{
    public class SingleSampleCommonOptionsParser : Option<SingleSampleCommonOptions>
    {
        private static readonly FileOption SampleBAlleleSites = FileOption.CreateRequired("vcf containing SNV b-allele sites in the sample (only sites with PASS in the filter column will be used)", "sample-b-allele-vcf");
        public static readonly FileOption PopulationBAlleleSites = FileOption.CreateRequired("vcf containing SNV b-allele sites in the population (only sites with PASS in the filter column will be used)", "population-b-allele-vcf");
        private static readonly ExclusiveFileOption BAlleleSites = ExclusiveFileOption.CreateRequired(SampleBAlleleSites, PopulationBAlleleSites);
        private static readonly FileOption PloidyBed = FileOption.Create(".bed file containing regions of known ploidy in the sample. Copy number calls matching the known ploidy in these regions will be considered non-variant", "ploidy-bed");
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