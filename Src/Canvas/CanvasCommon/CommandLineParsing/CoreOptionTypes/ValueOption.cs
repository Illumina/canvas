using System;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using System.Reflection;

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
                T parsedValue = ConvertInternal(value);
                return ParsingResult<T>.SuccessfulResult(parsedValue);
            }
            catch (Exception)
            {
                return ParsingResult<T>.FailedResult($"Error parsing {Info.Name} option: failed to convert {value} to {typeof(T)}");
            }
        }

        private static T ConvertInternal(string value)
        {
            if (typeof(T).GetTypeInfo().IsEnum)
                return (T)Enum.Parse(typeof(T), value, true);

            return (T)Convert.ChangeType(value, GetUnderlyingType(typeof(T)));
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

        public override bool Equals(object obj)
        {
            // If parameter is null return false.
            if (obj == null)
            {
                return false;
            }

            // If parameter cannot be cast to OptionInfo<T> return false.
            var p = obj as ValueOption<T>;
            if (p == null)
            {
                return false;
            }

            // Return true if the fields match:
            return Equals(p);
        }

        public bool Equals(ValueOption<T> p)
        {
            // If parameter is null return false:
            if (p == null)
            {
                return false;
            }

            return Info.Equals(p.Info) && Equals(_defaultValue, p._defaultValue);
        }

        public override int GetHashCode()
        {
            return Info.GetHashCode();
        }
    }
}