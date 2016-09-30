using System;
using CanvasCommon.CommandLineParsing.OptionProcessing;

namespace CanvasCommon.CommandLineParsing.CoreOptionTypes
{
    public class ValueOption<T> : Option<T>
    {
        private readonly T _defaultValue;
        public ValueOptionInfo<string> Info { get; }

        protected ValueOption(ValueOptionInfo<string> info, T defaultValue = default(T))
        {
            _defaultValue = defaultValue;
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

        public static ValueOption<T> CreateWithDefault(T defaultValue, string description, params string[] names)
        {
            return new ValueOption<T>(new ValueOptionInfo<string>(false, description, names), defaultValue);
        }

        public override ParsingResult<T> Parse(SuccessfulResultCollection value)
        {
            return Parse(value.Get(Info));
        }

        public virtual ParsingResult<T> Parse(string value)
        {
            if (value == null) return ParsingResult<T>.SuccessfulResult(_defaultValue);
            try
            {
                T parsedValue = (T)Convert.ChangeType(value, GetUnderlyingType(typeof(T)));
                return ParsingResult<T>.SuccessfulResult(parsedValue);
            }
            catch (Exception)
            {
                return ParsingResult<T>.FailedResult($"Error parsing {Info.Name} option: failed to convert {value} to {typeof(T)}");
            }
        }

        private static Type GetUnderlyingType(Type type)
        {
            Type underlyingType = Nullable.GetUnderlyingType(type);
            if (underlyingType == null)
                return type;
            return underlyingType;
        }

        public override OptionCollection<T> GetOptions()
        {
            return new OptionCollection<T> { Info };
        }
    }
}