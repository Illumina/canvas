using System.IO;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Isas.Framework.Checkpointing;
using Isas.Framework.FrameworkFactory;
using Isas.Framework.Logging;
using Isas.Framework.Settings;
using Isas.Framework.WorkManagement;
using static CanvasCommon.CommandLineParsing.CoreOptionTypes.OptionExtensions;

namespace Canvas.CommandLineParsing
{
    public interface IModeParser
    {
        string Name { get; }
        string Description { get; }
        IParsingResult<IModeLauncher> Parse(MainParser main, ILogger logger, ISettings settings, ICheckpointRunner checkpointRunner, IWorkManager workManager, IWorkDoer workDoer, string[] args, TextWriter standardWriter,
            TextWriter errorWriter);
        OptionCollection<IModeLauncher> GetModeOptions();
        void ShowHelp(WriteLine writer);
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

        public IParsingResult<IModeLauncher> Parse(MainParser main, ILogger logger, ISettings settings, ICheckpointRunner checkpointRunner, IWorkManager workManager, IWorkDoer workDoer, string[] args, TextWriter standardWriter, TextWriter errorWriter)
        {
            var results = main.GetParseResults(args);

            var parsingResult = checkpointRunner.RunCheckpoint("Validate input", () =>
            {
                if (!results.Validate(out IParsingResult<IModeLauncher> failedResult))
                    return ParsingResult<T>.FailedResult(failedResult.ErrorMessage);
                var successfulResults = new SuccessfulResultCollection(results);
                var commonOptions = successfulResults.Get(MainParser.CommonOptionsParser);
                return GetSerializedResult(successfulResults, commonOptions);
            });

            if (!parsingResult.Success)
            {
                ShowError(main, logger.Error, parsingResult.ErrorMessage);
                return ParsingResult<IModeLauncher>.FailedResult(parsingResult);
            }
            var runner = GetRunner(parsingResult.Result);
            return ParsingResult<IModeLauncher>.SuccessfulResult(new ModeLauncher(logger, settings, checkpointRunner, workManager, workDoer, runner, args, main.GetVersion(), Name));
        }

        private void ShowError(MainParser main, WriteLine errorWriter, string errorMessage)
        {
            errorWriter(errorMessage);
            errorWriter(" ");
            main.ShowHelp(errorWriter, this);
        }

        public void ShowHelp(WriteLine writer)
        {
            GetModeOptions().ShowHelp(writer);
        }

        public abstract IParsingResult<T> GetSerializedResult(SuccessfulResultCollection result, CommonOptions commonOptions);

        public abstract IModeRunner GetRunner(T result);

        public abstract OptionCollection<IModeLauncher> GetModeOptions();
    }
}