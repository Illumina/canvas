using Canvas.SmallPedigree;
using CanvasCommon.CommandLineParsing.OptionProcessing;

namespace Canvas.CommandLineParsing
{
    public class SmallPedigreeModeParser : ModeParser
    {
        private static readonly CommonOptionsParser CommonOptionsParser = new CommonOptionsParser();
        private static readonly SmallPedigreeOptionsParser SmallPedigreeOptionsParser = new SmallPedigreeOptionsParser();

        public SmallPedigreeModeParser(string name, string description) : base(name, description)
        {
        }

        public override ParsingResult<IModeRunner> Parse(SuccessfulResultCollection parseInput)
        {
            CommonOptions commonOptions = parseInput.Get(CommonOptionsParser);
            SmallPedigreeOptions smallPedigreeOptions = parseInput.Get(SmallPedigreeOptionsParser);
            return ParsingResult<IModeRunner>.SuccessfulResult(new SmallPedigreeWgsRunner(commonOptions, smallPedigreeOptions));
        }

        public override OptionCollection<IModeRunner> GetOptions()
        {
            return new OptionCollection<IModeRunner>
            {
                SmallPedigreeOptionsParser, CommonOptionsParser
            };
        }
    }
}