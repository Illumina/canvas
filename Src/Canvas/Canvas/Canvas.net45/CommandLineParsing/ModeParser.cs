using System;
using System.IO;
using System.Linq;
using CanvasCommon.CommandLineParsing.CoreOptionTypes;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Isas.Framework.FrameworkFactory;

namespace Canvas.CommandLineParsing
{
    public interface IModeParser : IOption
    {
        string Name { get; }
        string Description { get; }
        ParsingResult<IModeLauncher> Parse(MainParser main, FrameworkServices frameworkServices, string[] args, TextWriter standardWriter,
            TextWriter errorWriter);
        OptionCollection<IModeLauncher> GetModeOptions();
    }

    public abstract class ModeParser<T> : Option<IModeLauncher>, IModeParser
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
            bool useResult = false;
            CommonOptions commonOptions = null;
            ParsingResult<IModeLauncher> parsingResult = ParsingResult<IModeLauncher>.SuccessfulResult(new NullModeLauncher());
            var t = frameworkServices.Checkpointer.RunCheckpoint("Validate input", () =>
            {
                var options = GetModeOptions();
                options.Add(MainParser.BaseOptionsParser);
                options.Add(MainParser.CommonOptionsParser);
                var results = options.Parse(args.Skip(1));
                var baseOptions = results.Get(MainParser.BaseOptionsParser);
                if (!results.RemainingArgs.Any() && baseOptions.Success && main.HandleBaseOptions(baseOptions.Result, standardWriter, this))
                {
                    useResult = true;
                    return default(T);
                }

                if (!results.Validate(out ParsingResult<IModeLauncher> failedResult))
                {
                    errorWriter.WriteLine(failedResult.ErrorMessage);
                    errorWriter.WriteLine();
                    main.ShowHelp(errorWriter, this);
                    useResult = true;
                    parsingResult = ParsingResult<IModeLauncher>.FailedResult(failedResult.ErrorMessage);
                    return default(T);
                }
                var successfulResults = new SuccessfulResultCollection(results);
                commonOptions = successfulResults.Get(MainParser.CommonOptionsParser);
                var result = GetResult(successfulResults, commonOptions);
                if (!result.Success)
                {
                    parsingResult = ParsingResult<IModeLauncher>.FailedResult(result.ErrorMessage);
                    useResult = false;
                    return default(T);
                }
                return result.Result;
            });
            if (useResult) return parsingResult;
            var runnerResult = GetRunner(t);
            if (!runnerResult.Success) return ParsingResult<IModeLauncher>.FailedResult(runnerResult.ErrorMessage);
            return ParsingResult<IModeLauncher>.SuccessfulResult(new ModeLauncher(frameworkServices, commonOptions, runnerResult.Result, args, main.GetVersion(), Name));
        }

        public abstract ParsingResult<T> GetResult(SuccessfulResultCollection result, CommonOptions commonOptions);

        public abstract ParsingResult<IModeRunner> GetRunner(T result);

        public abstract OptionCollection<IModeLauncher> GetModeOptions();
        public override ParsingResult<IModeLauncher> Parse(SuccessfulResultCollection parseInput)
        {
            throw new NotSupportedException();
        }

        public override OptionCollection<IModeLauncher> GetOptions()
        {
            return GetModeOptions();
        }
    }
}