using System;
using Illumina.Common;

namespace CanvasCommon.CommandLineParsing.OptionProcessing
{
    public interface IParsingResult
    {
        bool Success { get; }
        string ErrorMessage { get; }
    }

    public class ParsingResult<T> : IParsingResult
    {
        public bool Success { get; }
        public string ErrorMessage { get; }
        public T Result { get; }

        private ParsingResult(bool success, T result, string errorMessage)
        {
            Result = result;
            Success = success;
            ErrorMessage = errorMessage;
        }
        public static ParsingResult<T> FailedResult(string errorMessage)
        {
            if (errorMessage.IsNullOrWhiteSpace())
            {
                throw new ArgumentException($"{nameof(errorMessage)} cannot be null or whitespace");
            }
            return new ParsingResult<T>(false, default(T), errorMessage);
        }
        public static ParsingResult<T> SuccessfulResult(T result)
        {
            return new ParsingResult<T>(true, result, "");
        }
    }
}