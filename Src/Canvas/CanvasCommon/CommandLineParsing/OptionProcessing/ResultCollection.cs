using System.Collections.Generic;
using System.Linq;
using CanvasCommon.CommandLineParsing.CoreOptionTypes;

namespace CanvasCommon.CommandLineParsing.OptionProcessing
{
    public interface IResultCollection
    {
        bool Validate(out IParsingResult result);
        IParsingResult<T> Get<T>(Option<T> option);
        IEnumerable<string> RemainingArgs { get; }
    }

    public class ResultCollection<TResult> : IResultCollection
    {
        private readonly IParsingResult<TResult> _failedResult;
        private readonly Dictionary<IOption, IParsingResult> _results;
        public IEnumerable<string> RemainingArgs { get; }

        public ResultCollection(Dictionary<IOption, IParsingResult> results, IEnumerable<string> remainingArgs = null)
        {
            remainingArgs = remainingArgs ?? Enumerable.Empty<string>();
            _results = results;
            RemainingArgs = remainingArgs;
        }

        public ResultCollection(IParsingResult<TResult> failedResult)
        {
            _failedResult = failedResult;
            RemainingArgs = Enumerable.Empty<string>();
        }

        public IParsingResult<T> Get<T>(Option<T> option)
        {
            if (_failedResult != null)
                return ParsingResult<T>.FailedResult(_failedResult.ErrorMessage);
            return (IParsingResult<T>)_results[option];
        }

        bool IResultCollection.Validate(out IParsingResult failedResult)
        {
            IParsingResult<TResult> failedResult2;
            bool b = Validate(out failedResult2);
            failedResult = failedResult2;
            return b;
        }

        public bool Validate(out IParsingResult<TResult> failedResult, bool allowUnparsedArguments = false)
        {
            if (_failedResult != null)
            {
                failedResult = _failedResult;
                return false;
            }
            if (!allowUnparsedArguments && RemainingArgs.Any())
            {
                failedResult = ParsingResult<TResult>.FailedResult($"Error: found unexpected arguments '{string.Join(" ", RemainingArgs)}'");
                return false;
            }
            failedResult = null;
            foreach (IParsingResult result in _results.Values)
            {
                if (result.Success) continue;
                failedResult = ParsingResult<TResult>.FailedResult(result.ErrorMessage);
                return false;
            }
            return true;
        }
    }
}