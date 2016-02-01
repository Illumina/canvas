using Canvas.CommandLineParsing.CoreOptionTypes;
using Canvas.CommandLineParsing.OptionProcessing;

namespace Canvas.CommandLineParsing
{
    public class GermlineWgsModeParser : ModeParser
    {
        private static readonly FileOption Bam = FileOption.CreateRequired("germline sample .bam file", "b", "bam");
        private static readonly CommonOptionsParser CommonOptionsParser = new CommonOptionsParser();

        public GermlineWgsModeParser(string name, string description) : base(name, description)
        {
        }

        public override ParsingResult<IModeRunner> Parse(SuccessfulResultCollection parseInput)
        {
            CommonOptions commonOptions = parseInput.Get(CommonOptionsParser);
            var bam = parseInput.Get(Bam);
            return ParsingResult<IModeRunner>.SuccesfulResult(new GermlineWgsRunner(commonOptions, bam));
        }

        public override OptionCollection<IModeRunner> GetOptions()
        {
            return new OptionCollection<IModeRunner> { Bam, CommonOptionsParser };
        }
    }
}