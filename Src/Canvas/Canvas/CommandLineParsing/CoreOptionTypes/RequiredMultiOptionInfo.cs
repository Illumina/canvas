using System.Collections.Generic;
using System.Linq;
using Canvas.CommandLineParsing.OptionProcessing;

namespace Canvas.CommandLineParsing.CoreOptionTypes
{
    public class RequiredMultiOptionInfo : ValueOptionInfo<List<string>>
    {
        public RequiredMultiOptionInfo(string description, params string[] names) : base(true, description, names)
        {
        }

        public RequiredMultiOptionInfo(IOptionInfo optionInfo) : base(true, optionInfo)
        { }


        public override ParsingResult<List<string>> Parse(SuccessfulResultCollection parseInput)
        {
            List<string> value = parseInput.Get(this);
            if (!value.Any())
                return ParsingResult<List<string>>.FailedResult($"Error: {Name} is a required option");
            return ParsingResult<List<string>>.SuccesfulResult(value);
        }
    }
}