using CanvasCommon.CommandLineParsing.OptionProcessing;

namespace CanvasCommon.CommandLineParsing.CoreOptionTypes
{
    public class StringOption : ValueOption<string>
    {
        private StringOption(ValueOptionInfo<string> option) : base(option)
        {
        }

        public static new StringOption CreateRequired(string description, params string[] names)
        {
            return new StringOption(new RequiredValueOptionInfo(description, names));
        }

        public static new StringOption Create(string description, params string[] names)
        {
            return new StringOption(new ValueOptionInfo<string>(false, description, names));
        }

        public override ParsingResult<string> Parse(string value)
        {
            return ParsingResult<string>.SuccessfulResult(value);
        }
    }
}