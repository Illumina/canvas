using System;
using Canvas.CommandLineParsing.OptionProcessing;

namespace Canvas.CommandLineParsing.CoreOptionTypes
{
    public class ValueOption<T> : Option<T>
    {
        public ValueOptionInfo<string> Info { get; }

        protected ValueOption(ValueOptionInfo<string> info)
        {
            Info = info;
        }

        public static ValueOption<T> CreateRequired(string description, params string[] names)
        {
            return new ValueOption<T>(new RequiredValueOptionInfo(description, names));
        }

        public static ValueOption<T> Create(string description, params string[] names)
        {
            return new ValueOption<T>(new ValueOptionInfo<string>(false, description, names));
        }

        public override ParsingResult<T> Parse(SuccessfulResultCollection value)
        {
            return Parse(value.Get(Info));
        }

        public virtual ParsingResult<T> Parse(string value)
        {
            if (value == null) return ParsingResult<T>.SuccesfulResult(default(T));
            try
            {
                T parsedValue = (T)Convert.ChangeType(value, typeof(T));
                return ParsingResult<T>.SuccesfulResult(parsedValue);
            }
            catch (Exception)
            {
                return ParsingResult<T>.FailedResult($"Error parsing {Info.Name} option: failed to convert {value} to {typeof(T)}");
            }
        }

        public override OptionCollection<T> GetOptions()
        {
            return new OptionCollection<T> { Info };
        }
    }
}