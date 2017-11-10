using System.Collections.Generic;
using System.Linq;
using CanvasCommon.CommandLineParsing.OptionProcessing;

namespace CanvasCommon.CommandLineParsing.CoreOptionTypes
{
    public class RequiredMultiOptionInfo : ValueOptionInfo<List<string>>
    {
        public RequiredMultiOptionInfo(string description, params string[] names) : base(true, description, names)
        {
        }

        public RequiredMultiOptionInfo(IOptionInfo optionInfo) : base(true, optionInfo)
        { }


        public override IParsingResult<List<string>> Parse(SuccessfulResultCollection parseInput)
        {
            List<string> value = parseInput.Get(this);
            if (!value.Any())
                return ParsingResult<List<string>>.FailedResult($"Error: {Name} is a required option");
            return ParsingResult<List<string>>.SuccessfulResult(value);
        }
    }
}