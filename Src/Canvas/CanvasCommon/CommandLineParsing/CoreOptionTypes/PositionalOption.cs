using System;
using System.Collections.Generic;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Illumina.Common;
using System.Linq;

namespace CanvasCommon.CommandLineParsing.CoreOptionTypes
{
    public class PositionalOption<T1, T2, T3, TOut> : Option<List<TOut>>
    {
        private ValueOption<T1> _option1;
        private ValueOption<T2> _option2;
        private ValueOption<T3> _option3;
        private Func<T1, T2, T3, IParsingResult<TOut>> _parse;
        private bool _required;

        public OptionInfo<List<List<string>>> Info { get; }
        public PositionalOption(Func<T1, T2, T3, IParsingResult<TOut>> parse, bool required, ValueOption<T1> option1, ValueOption<T2> option2, ValueOption<T3> option3, params string[] names)
        {
            _required = required;
            _parse = parse;
            _option1 = option1;
            _option2 = option2;
            _option3 = option3;
            Info = new PositionalOptionInfo(required, new IOptionInfo[] { option1.Info, option2.Info, option3.Info }, names);
        }

        public override IParsingResult<List<TOut>> Parse(SuccessfulResultCollection input)
        {
            var multipleValues = input.Get(Info);
            var outputs = new List<TOut>();
            foreach (var values in multipleValues)
            {
                var option1Result = GetOptionResult(0, values, _option1);
                var option2Result = GetOptionResult(1, values, _option2);
                var option3Result = GetOptionResult(2, values, _option3);
                if (!option1Result.Success) return ParsingResult<List<TOut>>.FailedResult(option1Result.ErrorMessage);
                if (!option2Result.Success) return ParsingResult<List<TOut>>.FailedResult(option2Result.ErrorMessage);
                if (!option3Result.Success) return ParsingResult<List<TOut>>.FailedResult(option3Result.ErrorMessage);
                var result = _parse(option1Result.Result, option2Result.Result, option3Result.Result);
                if (!result.Success) return ParsingResult<List<TOut>>.FailedResult(result.ErrorMessage);
                outputs.Add(result.Result);
            }
            if (_required && outputs.Empty())
                return ParsingResult<List<TOut>>.FailedResult($"{Info.Name} is a required option");
            return ParsingResult<List<TOut>>.SuccessfulResult(outputs);
        }

        private IParsingResult<T> GetOptionResult<T>(int valueIndex, List<string> values, ValueOption<T> option1)
        {
            string value = valueIndex < values.Count ? values[valueIndex] : null;

            if (option1.Info is RequiredValueOptionInfo && value == null)
                return ParsingResult<T>.FailedResult($"{option1.Info.Name} is a required positional argument for option {Info.Name}");
            return option1.Parse(value);
        }

        public override OptionCollection<List<TOut>> GetOptions()
        {
            return new OptionCollection<List<TOut>> { Info };
        }
    }

    public class PositionalOptionInfo : OptionInfo<List<List<string>>>
    {
        public int MaxNumValues { get; }
        private readonly bool _required;

        public PositionalOptionInfo(bool required, string description, int maxNumValues, params string[] names) : base(description, names)
        {
            MaxNumValues = maxNumValues;
            _required = required;
        }

        public PositionalOptionInfo(bool required, IOptionInfo[] optionsInfos, params string[] names) : this(false, GetDescription(required, optionsInfos), optionsInfos.Length, names)
        {
        }

        private static string GetDescription(bool required, params IOptionInfo[] optionsInfos)
        {
            var names = string.Join(" ", optionsInfos.Select(info =>
            {
                if (info is RequiredValueOptionInfo)
                    return info.Name;
                return $"[{info.Name}]";
            }));

            var descriptions = "";
            foreach (var option in optionsInfos)
            {
                if (descriptions.Length == 0) descriptions += " ";
                descriptions += $"{option.Name}: {option.Description}\n\n";
            }
            names += " Option can be specified multiple times.";
            if (required)
                names += " (required)";
            return names + "\n\n" + descriptions.Trim();
        }

        public override string Description
        {
            get
            {
                string description = base.Description;
                if (_required)
                {
                    description += " (required)";
                }
                return description;
            }
        }

        public override string GetPrototype()
        {
            return base.GetPrototype() + "={}";
        }
    }

    public sealed class PositionalOptionSetOption : Option
    {
        public delegate void PositionalAction(params string[] positionalArguments);
        PositionalAction action;

        public PositionalOptionSetOption(string prototype, string description, int count, PositionalAction action)
            : base(prototype, description, count)
        {
            if (action == null)
                throw new ArgumentNullException("action");
            this.action = action;
        }

        protected override void OnParseComplete(OptionContext c)
        {
            action(c.OptionValues.ToArray());
        }
    }

    public static class OptionSetExtensions
    {
        public static OptionSet Add(this OptionSet optionSet, string prototype, string description, int count, PositionalOptionSetOption.PositionalAction action)
        {
            return optionSet.Add(new PositionalOptionSetOption(prototype, description, count, action));
        }
    }
}