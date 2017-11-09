using System.Collections.Generic;
using System.Linq;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Illumina.Common;

namespace CanvasCommon.CommandLineParsing.CoreOptionTypes
{
    public class DictionaryOption : Option<Dictionary<string, string>>
    {
        private readonly MultiValueOption<string> _option;
        private readonly string _separator;

        private DictionaryOption(MultiValueOption<string> option, string separator)
        {
            _option = option;
            _separator = separator;
        }

        public static DictionaryOption Create(string description, params string[] names)
        {
            return new DictionaryOption(new MultiValueOption<string>(StringOption.Create(description, names)), ",");
        }

        public static DictionaryOption CreateRequired(string description, params string[] names)
        {
            return new DictionaryOption(new MultiValueOption<string>(StringOption.CreateRequired(description, names)), ",");
        }

        public override OptionCollection<Dictionary<string, string>> GetOptions()
        {
            return new OptionCollection<Dictionary<string, string>>
            {
                _option
            };
        }

        public override IParsingResult<Dictionary<string, string>> Parse(SuccessfulResultCollection parseInput)
        {
            List<string> inputs = parseInput.Get(_option);
            var result = new Dictionary<string, string>();
            foreach (var input in inputs)
            {
                string[] split = input.Split(_separator).Select(untrimmed => untrimmed.Trim()).ToArray();
                if (split.Length != 2)
                    return ParsingResult<Dictionary<string, string>>.FailedResult($"Error: expected format is {{KEY}}{_separator}{{value}}. Input was {input}");
                string key = split[0];
                string value = split[1];
                if (result.ContainsKey(key))
                {
                    return ParsingResult<Dictionary<string, string>>.FailedResult($"Error parsing {_option.Info.Name}: found duplicate key {key}");
                }
                result.Add(key, value);
            }
            return ParsingResult<Dictionary<string, string>>.SuccessfulResult(result);
        }
    }
}