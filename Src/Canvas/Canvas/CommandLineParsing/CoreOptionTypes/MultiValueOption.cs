using System.Collections.Generic;
using Canvas.CommandLineParsing.OptionProcessing;

namespace Canvas.CommandLineParsing.CoreOptionTypes
{
    public class MultiValueOption<T> : Option<List<T>>
    {
        private readonly ValueOption<T> _valueOption;
        public ValueOptionInfo<List<string>> Info { get; }

        public MultiValueOption(ValueOption<T> valueOption)
        {
            _valueOption = valueOption;
            if (valueOption.Info is RequiredValueOptionInfo)
            {
                Info = new RequiredMultiOptionInfo(valueOption.Info);
            }
            else
            {
                Info = new MultipleValueOptionInfo(false, valueOption.Info);
            }
        }

        public override ParsingResult<List<T>> Parse(SuccessfulResultCollection input)
        {
            var values = input.Get(Info);
            var results = new List<T>();
            foreach (var value in values)
            {
                var result = _valueOption.Parse(value);
                if (!result.Success)
                    return ParsingResult<List<T>>.FailedResult(result.ErrorMessage);
                results.Add(result.Result);
            }
            return ParsingResult<List<T>>.SuccesfulResult(results);
        }

        public override OptionCollection<List<T>> GetOptions()
        {
            return new OptionCollection<List<T>> { Info };
        }
    }
}