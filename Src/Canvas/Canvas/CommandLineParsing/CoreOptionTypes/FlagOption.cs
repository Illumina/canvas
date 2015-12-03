using Canvas.CommandLineParsing.OptionProcessing;

namespace Canvas.CommandLineParsing.CoreOptionTypes
{
    public class FlagOption : Option<bool>
    {
        private readonly OptionInfo<string> _info;

        public FlagOption(string description, params string[] names)
        {
            _info = new OptionInfo<string>(description, names);
        }

        public ParsingResult<bool> Parse(string value)
        {
            return ParsingResult<bool>.SuccesfulResult(value != null);
        }

        public override OptionCollection<bool> GetOptions()
        {
            return new OptionCollection<bool> { _info };
        }

        public override ParsingResult<bool> Parse(SuccessfulResultCollection parseInput)
        {
            return Parse(parseInput.Get(_info));
        }
    }
}