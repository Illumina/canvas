using CanvasCommon.CommandLineParsing.CoreOptionTypes;
using CanvasCommon.CommandLineParsing.OptionProcessing;

namespace Canvas.CommandLineParsing
{
    public class BaseOptionsParser : Option<BaseOptions>
    {
        private static readonly FlagOption Help = new FlagOption("show this message and exit", "h", "help");
        private static readonly FlagOption Version = new FlagOption("print version and exit", "v", "version");

        public override ParsingResult<BaseOptions> Parse(SuccessfulResultCollection parseInput)
        {
            var helpResult = parseInput.Get(Help);
            var versionResult = parseInput.Get(Version);
            return ParsingResult<BaseOptions>.SuccessfulResult(new BaseOptions(helpResult, versionResult));
        }

        public override OptionCollection<BaseOptions> GetOptions()
        {
            return new OptionCollection<BaseOptions>
            {
                Help,
                Version
            };
        }
    }
}