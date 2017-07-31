using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using CanvasCommon.CommandLineParsing.CoreOptionTypes;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Isas.Framework.FrameworkFactory;
using Isas.Framework.Settings;
using Isas.Framework.WorkManagement;

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

        public IParsingResult<IModeLauncher> Parse(FrameworkServices frameworkServices, string[] args, TextWriter standardWriter, TextWriter errorWriter)
        {
            var mode = _modeParsers[args[0]];
            return mode.Parse(this, frameworkServices, args, standardWriter, errorWriter);
        }

        public int Parse(string[] args, TextWriter standardWriter, TextWriter errorWriter)
        {
            var mode = _modeParsers[args[0]];
            var options = mode.GetModeOptions();
            options.Add(BaseOptionsParser);
            options.Add(CommonOptionsParser);
            var result = options.Parse(args.Skip(1));
            var baseOptions = result.Get(BaseOptionsParser);
            if (!result.RemainingArgs.Any() && baseOptions.Success && HandleBaseOptions(baseOptions.Result, standardWriter))
                return 0;

            if (!result.Validate(out IParsingResult<IModeLauncher> failedResult))
            {
                errorWriter.WriteLine(failedResult.ErrorMessage);
                errorWriter.WriteLine();
                ShowHelp(errorWriter);
                return -1;
            }

            throw new IlluminaException("This method should only be called when we fail to parse common options");
        }

        public int HandleMissingMode(IEnumerable<string> args, TextWriter standardWriter, TextWriter errorWriter)
        {
            string error = "Error: no mode specified";
            var baseOptionsResult = BaseOptionsParser.Parse(args);
            if (baseOptionsResult.Success)
            {
                if (HandleBaseOptions(baseOptionsResult.Result, standardWriter))
                    return 0;
            }
            else
            {
                error = baseOptionsResult.ErrorMessage;
            }
            errorWriter.WriteLine(error);
            errorWriter.WriteLine();
            ShowHelp(errorWriter);
            return -1;
        }

        public bool HandleBaseOptions(BaseOptions baseOptions, TextWriter writer, IModeParser specifiedMode = null)
        {
            if (baseOptions.ShowHelp)
            {
                ShowHelp(writer, specifiedMode);
                return true;
            }

            if (baseOptions.ShowVersion)
            {
                ShowVersion(writer);
                return true;
            }
            return false;
        }

        public void ShowHelp(TextWriter writer, IModeParser specifiedMode = null)
        {
            writer.WriteLine($"Canvas {GetVersion()} {GetCopyright()}");
            writer.WriteLine();
            string modeName = specifiedMode?.Name ?? "[MODE]";
            if (specifiedMode != null)
            {
                writer.WriteLine($"{specifiedMode.Name} - {specifiedMode.Description}");
                writer.WriteLine();
            }
            writer.WriteLine($"Usage: Canvas.exe {modeName} [OPTIONS]+");
            writer.WriteLine();
            if (specifiedMode == null)
            {
                writer.WriteLine("Available modes:");
                foreach (var mode in _modeParsers.Values)
                {
                    writer.WriteLine($"\t{mode.Name} - {mode.Description}");
                }
                writer.WriteLine();
            }
            writer.WriteLine("Options:");
            if (specifiedMode != null)
            {
                specifiedMode.ShowHelp(writer);
            }
            BaseOptionsParser.ShowHelp(writer);
        }

        private void ShowVersion(TextWriter standardWriter)
        {
            standardWriter.WriteLine(GetVersion());
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
                return Parse(args, standardOutput, standardError);
            }
            var commonOptions = commonOptionsResult.Result;
            IDirectoryLocation outFolder = commonOptions.OutputDirectory;
            var log = outFolder.GetFileLocation("CanvasLog.txt");
            var error = outFolder.GetFileLocation("CanvasError.txt");
            IsasConfiguration config = IsasConfiguration.GetConfiguration();
            IDirectoryLocation genomeRoot = commonOptions.WholeGenomeFasta?.Parent?.Parent?.Parent?.Parent?.Parent;
            int returnValue = 0;
            var settings = new SettingsProcessor();
            var maximumMemoryGB = settings.GetSetting(WorkManagerFactory.SampleSheetSettings.MaximumMemoryGB);
            var maximumHoursPerProcess = settings.GetSetting(WorkManagerFactory.SampleSheetSettings.MaximumHoursPerProcess);
            var maximumThreadCount = settings.GetSetting(WorkManagerFactory.SampleSheetSettings.MaximumThreadCount);
            IsasFrameworkFactory.RunWithIsasFramework(outFolder, log, error, commonOptions.StartCheckpoint, commonOptions.StopCheckpoint, maximumThreadCount,
                maximumMemoryGB, maximumHoursPerProcess, true, false, true, genomeRoot,
                frameworkServices =>
                {
                    var result = Parse(frameworkServices, args, standardOutput, standardError);
                    if (!result.Success)
                        returnValue = -1;
                    returnValue = result.Result.Launch();
                });
            return returnValue;
        }
    }
}