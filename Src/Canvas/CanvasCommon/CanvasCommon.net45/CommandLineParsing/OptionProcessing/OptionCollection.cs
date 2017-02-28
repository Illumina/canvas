using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using CanvasCommon.CommandLineParsing.CoreOptionTypes;
using NDesk.Options;

namespace CanvasCommon.CommandLineParsing.OptionProcessing
{
    public class OptionData
    {
        public List<string> Data { get; } = new List<string>();
    }
    public interface IOptionCollection
    {
        IResultCollection ParseInternal(IEnumerable<string> args);
        IResultCollection ParseInternal(Dictionary<IOption, OptionData> optionResults, IEnumerable<string> remainingArgs = null);
        IEnumerable<IOptionInfo> GetLeafs();
    }

    public class OptionCollection<TResult> : IOptionCollection, IEnumerable<KeyValuePair<IOption, IOptionCollection>>
    {
        private IResultCollection ParseInternal(IEnumerable<string> args)
        {
            return Parse(args);
        }

        private IResultCollection ParseInternal(Dictionary<IOption, OptionData> optionResults, IEnumerable<string> remainingArgs)
        {
            return ParseGeneric(optionResults, remainingArgs);
        }

        private readonly Dictionary<IOption, IOptionCollection> _optionCollections = new Dictionary<IOption, IOptionCollection>();
        public void Add<T>(Option<T> option)
        {
            IOptionInfo info = option as IOptionInfo;
            if (info == null)
            {
                _optionCollections.Add(option, option.GetOptions());
            }
            else
            {
                _optionCollections.Add(option, new OptionCollection<T>());
            }
        }

        public ResultCollection<TResult> Parse(IEnumerable<string> args)
        {
            Dictionary<IOption, OptionData> optionResults;
            OptionSet options = GetOptionSet(out optionResults);
            try
            {
                IEnumerable<string> remainingArgs = options.Parse(args);
                return ParseGeneric(optionResults, remainingArgs);
            }
            catch (OptionException e)
            {
                return new ResultCollection<TResult>(ParsingResult<TResult>.FailedResult(e.Message));
            }
        }

        private OptionSet GetOptionSet(out Dictionary<IOption, OptionData> optionResults)
        {
            optionResults = new Dictionary<IOption, OptionData>();
            OptionSet options = new OptionSet();
            foreach (var option in GetLeafs())
            {
                OptionData result = AddOption(option, options);
                optionResults.Add(option, result);
            }
            return options;
        }

        private OptionSet GetOptionSet()
        {
            Dictionary<IOption, OptionData> optionsResults;
            return GetOptionSet(out optionsResults);
        }

        private ResultCollection<TResult> ParseGeneric(Dictionary<IOption, OptionData> optionResults, IEnumerable<string> remainingArgs)
        {
            var results = new Dictionary<IOption, IParsingResult>();
            foreach (var option in _optionCollections)
            {
                if (optionResults.ContainsKey(option.Key))
                {
                    var baseOption = option.Key;
                    var optionData = optionResults[option.Key];
                    if (baseOption is OptionInfo<string>)
                        results.Add(option.Key, GetParseResult((OptionInfo<string>)baseOption, optionData));
                    else if (baseOption is OptionInfo<List<string>>)
                        results.Add(option.Key, GetParseResult((OptionInfo<List<string>>)baseOption, optionData));
                    continue;
                }
                IResultCollection resultCollection = option.Value.ParseInternal(optionResults);
                IParsingResult failedResult;
                if (!resultCollection.Validate(out failedResult))
                {
                    results.Add(option.Key, failedResult);
                }
                else
                {
                    results.Add(option.Key, option.Key.Parse(new SuccessfulResultCollection(resultCollection)));
                }
            }

            return new ResultCollection<TResult>(results, remainingArgs);
        }

        private static IParsingResult GetParseResult(OptionInfo<string> optionInfo, OptionData optionData)
        {
            ParsingResult<string> result = ParsingResult<string>.SuccessfulResult(optionData.Data.FirstOrDefault());
            if (optionData.Data.Count > 1)
                result = ParsingResult<string>.FailedResult($"Error: {optionInfo.Name} can only be specified once");
            return optionInfo.Parse(new SuccessfulResultCollection(optionInfo, result));
        }

        private static IParsingResult GetParseResult(OptionInfo<List<string>> multiOptionInfo, OptionData optionData)
        {
            return multiOptionInfo.Parse(new SuccessfulResultCollection(multiOptionInfo, ParsingResult<List<string>>.SuccessfulResult(optionData.Data)));
        }

        private OptionData AddOption(IOptionInfo info, OptionSet set)
        {
            OptionData s = new OptionData();
            set.Add(info.GetPrototype(), info.Description, v => s.Data.Add(v));
            return s;
        }

        private IEnumerable<IOptionInfo> GetLeafs()
        {
            var leafs = new List<IOptionInfo>();
            foreach (var option in _optionCollections)
            {

                IOptionInfo info = option.Key as IOptionInfo;
                if (option.Key is IOptionInfo)
                    leafs.Add(info);
                else
                    leafs.AddRange(option.Value.GetLeafs());
            }
            return leafs;
        }

        public void ShowHelp(TextWriter writer)
        {
            GetOptionSet().WriteOptionDescriptions(writer);
        }

        IResultCollection IOptionCollection.ParseInternal(IEnumerable<string> args)
        {
            return ParseInternal(args);
        }

        IResultCollection IOptionCollection.ParseInternal(Dictionary<IOption, OptionData> optionResults, IEnumerable<string> remainingArgs)
        {
            return ParseInternal(optionResults, remainingArgs);
        }

        IEnumerable<IOptionInfo> IOptionCollection.GetLeafs()
        {
            return GetLeafs();
        }

        public IEnumerator<KeyValuePair<IOption, IOptionCollection>> GetEnumerator()
        {
            return _optionCollections.GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return _optionCollections.GetEnumerator();
        }
    }
}