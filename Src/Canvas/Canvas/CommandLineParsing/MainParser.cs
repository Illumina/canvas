using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Canvas.CommandLineParsing.CoreOptionTypes;
using Canvas.CommandLineParsing.OptionProcessing;
using ILMNcommon.Common;

namespace Canvas.CommandLineParsing
{
    public class MainParser
    {
        private static readonly BaseOptionsParser BaseOptionsParser = new BaseOptionsParser();
        private readonly Dictionary<string, ModeParser> _modeParsers;

        public MainParser(params ModeParser[] modeParsers)
        {
            _modeParsers = modeParsers.ToDictionary(modeParser => modeParser.Name, modeParser => modeParser, StringComparer.OrdinalIgnoreCase);
        }

        public ParsingResult<IModeLauncher> Parse(string[] args, TextWriter standardWriter, TextWriter errorWriter)
        {
            if (args.Empty() || !_modeParsers.ContainsKey(args[0]))
            {
                return GetMissingModeResult(args, standardWriter, errorWriter);
            }
            var mode = _modeParsers[args[0]];
            var options = new OptionCollection<IModeLauncher> { BaseOptionsParser, mode };
            var result = options.Parse(args.Skip(1));
            ParsingResult<IModeLauncher> failedResult;
            if (!result.RemainingArgs.Any() && HandleBaseOptions(result.Get(BaseOptionsParser).Result, standardWriter, mode))
                return ParsingResult<IModeLauncher>.SuccesfulResult(new NullModeLauncher());

            if (!result.Validate(out failedResult))
            {
                errorWriter.WriteLine(failedResult.ErrorMessage);
                errorWriter.WriteLine();
                ShowHelp(errorWriter, mode);
                return failedResult;
            }
            var runner = result.Get(mode).Result;
            return ParsingResult<IModeLauncher>.SuccesfulResult(new ModeLauncher(runner, args, GetVersion(), mode.Name));
        }

        private ParsingResult<IModeLauncher> GetMissingModeResult(IEnumerable<string> args, TextWriter standardWriter, TextWriter errorWriter)
        {
            string error = "Error: no mode specified";
            var baseOptionsResult = BaseOptionsParser.Parse(args);
            if (baseOptionsResult.Success)
            {
                if (HandleBaseOptions(baseOptionsResult.Result, standardWriter))
                    return ParsingResult<IModeLauncher>.SuccesfulResult(new NullModeLauncher());
            }
            else
            {
                error = baseOptionsResult.ErrorMessage;
            }

            errorWriter.WriteLine(error);
            errorWriter.WriteLine();
            ShowHelp(errorWriter);
            return ParsingResult<IModeLauncher>.FailedResult(error);
        }

        private bool HandleBaseOptions(BaseOptions baseOptions, TextWriter writer, ModeParser specifiedMode = null)
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

        private void ShowHelp(TextWriter writer, ModeParser specifiedMode = null)
        {
            writer.WriteLine($"Canvas {GetVersion()} {Illumina.Shared.Version.VersionInfo.AssemblyCopyright}");
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

        private string GetVersion()
        {
            return typeof(MainParser).Assembly.GetName().Version.ToString();
        }
    }
}