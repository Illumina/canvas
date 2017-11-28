using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using CanvasCommon.CommandLineParsing.CoreOptionTypes;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Isas.Framework.Checkpointing;
using Isas.Framework.FrameworkFactory;
using Isas.Framework.Logging;
using Isas.Framework.Settings;
using Isas.Framework.WorkManagement;
using static CanvasCommon.CommandLineParsing.CoreOptionTypes.OptionExtensions;

namespace Canvas.CommandLineParsing
{
    public class MainParser
    {
        private readonly string _version;
        private readonly string _copyright;
        public static readonly BaseOptionsParser BaseOptionsParser = new BaseOptionsParser();
        public static readonly CommonOptionsParser CommonOptionsParser = new CommonOptionsParser();
        private readonly Dictionary<string, IModeParser> _modeParsers;

        public MainParser(string version, string copyright, params IModeParser[] modeParsers)
        {
            _version = version;
            _copyright = copyright;
            _modeParsers = modeParsers.ToDictionary(modeParser => modeParser.Name, modeParser => modeParser, StringComparer.OrdinalIgnoreCase);
        }

        public bool IsMissingMode(string[] args)
        {
            return args.Empty() || !_modeParsers.ContainsKey(args[0]);
        }

        public IParsingResult<CommonOptions> ParseCommonOptions(string[] args)
        {
            var options = new OptionCollection<IModeLauncher> { CommonOptionsParser };
            var result = options.Parse(args.Skip(1));
            return result.Get(CommonOptionsParser);
        }

        private IParsingResult<IModeLauncher> Parse(ILogger logger, ISettings settings, ICheckpointRunner checkpointRunner, IWorkManager workManager, IWorkDoer workDoer, string[] args, TextWriter standardWriter, TextWriter errorWriter)
        {
            var mode = _modeParsers[args[0]];
            return mode.Parse(this, logger, settings, checkpointRunner, workManager, workDoer, args, standardWriter, errorWriter);
        }

        private int Parse(string[] args, WriteLine standardWriter, WriteLine errorWriter)
        {
            var mode = GetMode(args);
            if (AnyBaseOptions(args, out BaseOptions baseOptions))
            {
                HandleBaseOptions(baseOptions, standardWriter, mode);
                return 0;
            }

            var result = GetParseResults(args);
            if (!result.Validate(out IParsingResult<IModeLauncher> failedResult))
            {
                errorWriter(failedResult.ErrorMessage);
                errorWriter(" ");
                ShowHelp(errorWriter, mode);
                return -1;
            }

            throw new IlluminaException("This method should only be called when we fail to parse common options");
        }

        private int HandleMissingMode(IEnumerable<string> args, TextWriter standardWriter, TextWriter errorWriter)
        {
            var baseOptionsResult = BaseOptionsParser.Parse(args);
            if (baseOptionsResult.Success && AnyBaseOptions(baseOptionsResult.Result))
            {
                HandleBaseOptions(baseOptionsResult.Result, standardWriter.WriteLine);
                return 0;
            }

            if (baseOptionsResult.Success)
            {
                errorWriter.WriteLine("Error: no mode specified");
            }
            else
            {
                errorWriter.WriteLine(baseOptionsResult.ErrorMessage);
            }
            errorWriter.WriteLine();
            ShowHelp(errorWriter.WriteLine);
            return -1;
        }

        private IModeParser GetMode(string[] args) => _modeParsers[args[0]];

        private bool AnyBaseOptions(string[] args, out BaseOptions baseOptionsResult)
        {
            baseOptionsResult = null;
            var result = GetParseResults(args);
            var baseOptions = result.Get(BaseOptionsParser);
            if (!result.RemainingArgs.Any() && baseOptions.Success && AnyBaseOptions(baseOptions.Result))
            {
                baseOptionsResult = baseOptions.Result;
                return true;
            }
            return false;
        }

        public ResultCollection<IModeLauncher> GetParseResults(string[] args)
        {
            var mode = GetMode(args);
            var options = mode.GetModeOptions();
            options.Add(BaseOptionsParser);
            options.Add(CommonOptionsParser);
            var result = options.Parse(args.Skip(1));
            return result;
        }

        private static bool AnyBaseOptions(BaseOptions baseOptions)
        {
            return baseOptions.ShowHelp || baseOptions.ShowVersion;
        }

        private void HandleBaseOptions(BaseOptions baseOptions, WriteLine writer, IModeParser specifiedMode = null)
        {
            if (baseOptions.ShowHelp)
            {
                ShowHelp(writer, specifiedMode);
            }

            if (baseOptions.ShowVersion)
            {
                ShowVersion(writer);
            }
        }

        public void ShowHelp(WriteLine writeLine, IModeParser specifiedMode = null)
        {
            writeLine($"Canvas {GetVersion()} {GetCopyright()}");
            writeLine(" ");
            string modeName = specifiedMode?.Name ?? "[MODE]";
            if (specifiedMode != null)
            {
                writeLine($"{specifiedMode.Name} - {specifiedMode.Description}");
                writeLine(" ");
            }
            writeLine($"Usage: Canvas.exe {modeName} [OPTIONS]+");
            writeLine(" ");
            if (specifiedMode == null)
            {
                writeLine("Available modes:");
                foreach (var mode in _modeParsers.Values)
                {
                    writeLine($"\t{mode.Name} - {mode.Description}");
                }
                writeLine(" ");
                writeLine("Options:");
                BaseOptionsParser.ShowHelp(writeLine);
            }
            else
            {
                writeLine("Mode-specific options:");
                specifiedMode.ShowHelp(writeLine);
                writeLine(" ");
                writeLine("Common options:");
                CommonOptionsParser.ShowHelp(writeLine);
                writeLine(" ");
                writeLine("Other options:");
                BaseOptionsParser.ShowHelp(writeLine);
            }
        }

        private void ShowVersion(WriteLine standardWriter)
        {
            standardWriter(GetVersion());
        }

        public string GetVersion()
        {
            return _version;
        }

        private string GetCopyright()
        {
            return _copyright;
        }

        public int Run(string[] args, TextWriter standardOutput, TextWriter standardError)
        {
            if (IsMissingMode(args))
            {
                return HandleMissingMode(args, standardOutput, standardError);
            }

            var commonOptionsResult = ParseCommonOptions(args);
            if (!commonOptionsResult.Success)
            {
                return Parse(args, standardOutput.WriteLine, standardError.WriteLine);
            }
            var commonOptions = commonOptionsResult.Result;

            if (AnyBaseOptions(args, out BaseOptions baseOptions))
            {
                HandleBaseOptions(baseOptions, standardOutput.WriteLine, GetMode(args));
                return 0;
            }
            IDirectoryLocation outFolder = commonOptions.OutputDirectory;
            var log = outFolder.GetFileLocation("CanvasLog.txt");
            var error = outFolder.GetFileLocation("CanvasError.txt");
            IDirectoryLocation genomeRoot = commonOptions.WholeGenomeFasta?.Parent?.Parent?.Parent?.Parent?.Parent;
            int returnValue = 0;

            IsasFrameworkFactory.RunWithLogger(log, error, logger =>
            {
                ISettings settings = new ProgrammaticSettingsBuilder()
                    .WithSetting(IsasConfigurationSettings.GenomeRepositoryPath, genomeRoot)
                    .ToSettings();
                settings = new NestedSettingsProcessor(settings, IsasConfigurationSettings.GetConfigSettings());
                IsasFrameworkFactory.RunWithWorkDoer(logger, settings, outFolder, workDoer =>
                {
                    IsasFrameworkFactory.RunWithCheckpointer(logger, outFolder, settings, commonOptions.StartCheckpoint,
                        commonOptions.StopCheckpoint, checkpointer =>
                    {
                        var workManager = WorkManagerFactory.GetWorkManager(workDoer, logger, outFolder, settings);
                        var result = Parse(logger, settings, checkpointer, workManager, workDoer, args, standardOutput, standardError);
                        returnValue = result.Success ? result.Result.Launch() : -1;
                    });
                });
            });

            return returnValue;
        }
    }
}