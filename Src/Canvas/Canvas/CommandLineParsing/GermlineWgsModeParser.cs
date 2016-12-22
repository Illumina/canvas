using CanvasCommon.CommandLineParsing.CoreOptionTypes;
using CanvasCommon.CommandLineParsing.OptionProcessing;

namespace Canvas.CommandLineParsing
{
    public class GermlineWgsModeParser : ModeParser
    {
        private static readonly FileOption Bam = FileOption.CreateRequired("germline sample .bam file", "b", "bam");
        private static readonly CommonOptionsParser CommonOptionsParser = new CommonOptionsParser();
        private static readonly SingleSampleCommonOptionsParser SingleSampleCommonOptionsParser = new SingleSampleCommonOptionsParser();

        public GermlineWgsModeParser(string name, string description) : base(name, description)
        {
        }

        public override ParsingResult<IModeRunner> Parse(SuccessfulResultCollection parseInput)
        {
            CommonOptions commonOptions = parseInput.Get(CommonOptionsParser);
            SingleSampleCommonOptions singleSampleCommonOptions = parseInput.Get(SingleSampleCommonOptionsParser);
            var bam = parseInput.Get(Bam);
            return ParsingResult<IModeRunner>.SuccessfulResult(new GermlineWgsRunner(commonOptions, singleSampleCommonOptions, bam));
        }

        public override OptionCollection<IModeRunner> GetOptions()
        {
            return new OptionCollection<IModeRunner> { Bam, CommonOptionsParser, SingleSampleCommonOptionsParser  };
        }
    }
}