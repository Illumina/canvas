using System;
using System.Linq;
using CanvasCommon.CommandLineParsing.OptionProcessing;

namespace CanvasCommon.CommandLineParsing.CoreOptionTypes
{
    public class SinglePositionalOption<T1, T2, T3, TOut> : Option<TOut>
    {
        public SinglePositionalOption(Func<T1, T2, T3, IParsingResult<TOut>> parse, bool required, string overallDescription, ValueOption<T1> option1, ValueOption<T2> option2, ValueOption<T3> option3, params string[] names)
        {
            PositionalOption = new PositionalOption<T1, T2, T3, TOut>(parse, false, required, overallDescription, option1, option2, option3, names);
        }

        private PositionalOption<T1, T2, T3, TOut> PositionalOption { get; }
        public override OptionCollection<TOut> GetOptions()
        {
            return new OptionCollection<TOut>
            {
                PositionalOption
            };
        }

        public override IParsingResult<TOut> Parse(SuccessfulResultCollection parseInput)
        {
            var multipleValues = parseInput.Get(PositionalOption);
            if (multipleValues.Count > 1)
                return ParsingResult<TOut>.FailedResult($"{PositionalOption.Info.Name} can be specified only once");
            return ParsingResult<TOut>.SuccessfulResult(multipleValues.SingleOrDefault());
        }
    }
}