using System.Collections.Generic;
using System.IO;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Isas.Framework.Logging;
using static CanvasCommon.CommandLineParsing.CoreOptionTypes.OptionExtensions;

namespace CanvasCommon.CommandLineParsing.CoreOptionTypes
{
    public interface IOption
    {
        IOptionCollection GetOptions();
        IParsingResult Parse(SuccessfulResultCollection parseInput);
        void ShowHelp(WriteLine writer);
    }

    public abstract class Option<TParseOutput> : IOption
    {
        public abstract OptionCollection<TParseOutput> GetOptions();
        public abstract IParsingResult<TParseOutput> Parse(SuccessfulResultCollection parseInput);
        public void ShowHelp(WriteLine writer)
        {
            OptionExtensions.ShowHelp(this, writer);
        }

        IOptionCollection IOption.GetOptions()
        {
            return GetOptions();
        }

        IParsingResult IOption.Parse(SuccessfulResultCollection parseInput)
        {
            return Parse(parseInput);
        }

    }

    public static class OptionExtensions
    {
        public delegate void WriteLine(string message);
        public static IParsingResult<T> Parse<T>(this Option<T> option, IEnumerable<string> args, bool allowUnparsedArguments = false)
        {
            OptionCollection<T> options = new OptionCollection<T>();
            options.Add(option);
            ResultCollection<T> result = options.Parse(args);
            IParsingResult<T> failedResult;
            if (!result.Validate(out failedResult, allowUnparsedArguments))
            {
                return failedResult;
            }
            return result.Get(option);
        }

        public static void ShowHelp<T>(this Option<T> option, WriteLine writer)
        {
            OptionCollection<T> options = new OptionCollection<T>();
            options.Add(option);
            options.ShowHelp(writer);
        }
    }
}