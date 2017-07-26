using System;
using System.IO;
using System.Linq;
using CanvasCommon.CommandLineParsing.CoreOptionTypes;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Isas.Framework.FrameworkFactory;

namespace Canvas.CommandLineParsing
{
    public interface IModeParser
    {
        string Name { get; }
        string Description { get; }
        ParsingResult<IModeLauncher> Parse(MainParser main, FrameworkServices frameworkServices, string[] args, TextWriter standardWriter,
            TextWriter errorWriter);
        OptionCollection<IModeLauncher> GetModeOptions();
        void ShowHelp(TextWriter writer);
    }

    public abstract class ModeParser<T> : IModeParser
    {
        public string Name { get; }
        public string Description { get; }

        protected ModeParser(string name, string description)
        {
            Name = name;
            Description = description;
        }

        public ParsingResult<IModeLauncher> Parse(MainParser main, FrameworkServices frameworkServices, string[] args, TextWriter standardWriter, TextWriter errorWriter)
        {
            var options = GetModeOptions();
            options.Add(MainParser.BaseOptionsParser);
            options.Add(MainParser.CommonOptionsParser);
            var results = options.Parse(args.Skip(1));
            var baseOptions = results.Get(MainParser.BaseOptionsParser);
            if (!results.RemainingArgs.Any() && baseOptions.Success && main.HandleBaseOptions(baseOptions.Result, standardWriter, this))
            {
                return ParsingResult<IModeLauncher>.SuccessfulResult(new NullModeLauncher());
            }

            var parsingResult = frameworkServices.Checkpointer.RunCheckpoint("Validate input", () =>
            {
                if (!results.Validate(out ParsingResult<IModeLauncher> failedResult))
                {
                    errorWriter.WriteLine(failedResult.ErrorMessage);
                    errorWriter.WriteLine();
                    main.ShowHelp(errorWriter, this);
                    return ParsingResult<T>.FailedResult(failedResult.ErrorMessage);
                }
                var successfulResults = new SuccessfulResultCollection(results);
                var commonOptions = successfulResults.Get(MainParser.CommonOptionsParser);
                return GetSerializedResult(successfulResults, commonOptions);
            });

            if (!parsingResult.Success) return ParsingResult<IModeLauncher>.FailedResult(parsingResult);
            var runner = GetRunner(parsingResult.Result);
            return ParsingResult<IModeLauncher>.SuccessfulResult(new ModeLauncher(frameworkServices, runner, args, main.GetVersion(), Name));
        }

        public void ShowHelp(TextWriter writer)
        {
            GetModeOptions().ShowHelp(writer);
        }

        public abstract ParsingResult<T> GetSerializedResult(SuccessfulResultCollection result, CommonOptions commonOptions);

        public abstract IModeRunner GetRunner(T result);

        public abstract OptionCollection<IModeLauncher> GetModeOptions();
    }
}