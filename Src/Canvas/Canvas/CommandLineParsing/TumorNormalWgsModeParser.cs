using Canvas.CommandLineParsing.CoreOptionTypes;
using Canvas.CommandLineParsing.OptionProcessing;
using Isas.Shared;

namespace Canvas.CommandLineParsing
{
    public class TumorNormalWgsModeParser : ModeParser
    {
        private static readonly CommonOptionsParser CommonOptionsParser = new CommonOptionsParser();
        private static readonly TumorNormalOptionsParser TumorNormalOptionsParser = new TumorNormalOptionsParser();

        public TumorNormalWgsModeParser(string name, string description) : base(name, description)
        {
        }

        public override ParsingResult<IModeRunner> Parse(SuccessfulResultCollection parseInput)
        {
            CommonOptions commonOptions = parseInput.Get(CommonOptionsParser);
            TumorNormalOptions tumorNormalOptions = parseInput.Get(TumorNormalOptionsParser);
            return ParsingResult<IModeRunner>.SuccesfulResult(new TumorNormalWgsRunner(commonOptions, tumorNormalOptions));
        }

        public override OptionCollection<IModeRunner> GetOptions()
        {
            return new OptionCollection<IModeRunner>
            {
                TumorNormalOptionsParser, CommonOptionsParser
            };
        }
    }

    internal class TumorNormalOptionsParser : Option<TumorNormalOptions>
    {
        private static readonly FileOption TumorBam = FileOption.CreateRequired("tumor sample .bam file", "b", "bam");
        private static readonly FileOption SomaticVcf = FileOption.Create(".vcf file of somatic small variants found in the tumor sample. " +
                                                                          "Used as a fallback for estimating the reported tumor sample purity. " +
                                                                          "This option has no effect on called CNVs", "somatic-vcf");

        public override OptionCollection<TumorNormalOptions> GetOptions()
        {
            return new OptionCollection<TumorNormalOptions>
            {
               TumorBam, SomaticVcf
            };
        }

        public override ParsingResult<TumorNormalOptions> Parse(SuccessfulResultCollection parseInput)
        {
            var tumorBam = parseInput.Get(TumorBam);
            var somaticVcf = parseInput.Get(SomaticVcf);

            return ParsingResult<TumorNormalOptions>.SuccesfulResult(
                new TumorNormalOptions(tumorBam, somaticVcf));
        }
    }

    public class TumorNormalOptions
    {
        public IFileLocation TumorBam { get; }
        public IFileLocation SomaticVcf { get; }

        public TumorNormalOptions(IFileLocation tumorBam, IFileLocation somaticVcf)
        {
            TumorBam = tumorBam;
            SomaticVcf = somaticVcf;
        }
    }
}