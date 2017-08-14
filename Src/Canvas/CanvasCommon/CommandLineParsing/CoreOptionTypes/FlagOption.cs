using CanvasCommon.CommandLineParsing.OptionProcessing;

namespace CanvasCommon.CommandLineParsing.CoreOptionTypes
{
    public class FlagOption : Option<bool>
    {
        private readonly OptionInfo<string> _info;

        public FlagOption(string description, params string[] names)
        {
            _info = new OptionInfo<string>(description, names);
        }

        public IParsingResult<bool> Parse(string value)
        {
            return ParsingResult<bool>.SuccessfulResult(value != null);
        }

        public override OptionCollection<bool> GetOptions()
        {
            return new OptionCollection<bool> { _info };
        }

        public override IParsingResult<bool> Parse(SuccessfulResultCollection parseInput)
        {
            return Parse(parseInput.Get(_info));
        }
    }
}