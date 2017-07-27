using System;
using Illumina.Common;

namespace CanvasCommon.CommandLineParsing.OptionProcessing
{
    public interface IParsingResult
    {
        bool Success { get; }
        string ErrorMessage { get; }
    }

    public interface IParsingResult<out T> : IParsingResult
    {
        T Result { get; }
    }

    public static class ParsingResult<T>
    {
        public static IParsingResult<T> FailedResult<U>(IParsingResult<U> failedResult)
        {
            if (failedResult.Success)
            {
                throw new ArgumentException($"{nameof(failedResult)} must not be successful");
            }
            return FailedResult(failedResult.ErrorMessage);
        }

        public static IParsingResult<T> FailedResult(string errorMessage)
        {
            if (string.IsNullOrWhiteSpace(errorMessage))
            {
                throw new ArgumentException($"{nameof(errorMessage)} cannot be null or whitespace");
            }
            return new FailedParsingResult<T>(errorMessage);
        }
        public static IParsingResult<T> SuccessfulResult(T result)
        {
            return new SuccessfulParsingResult<T>(result);
        }
    }

    public class FailedParsingResult<T> : IParsingResult<T>
    {
        public FailedParsingResult(string errorMessage)
        {
            ErrorMessage = errorMessage;
        }

        public bool Success => false;
        public string ErrorMessage { get; }
        public T Result => default(T);
    }

    public class SuccessfulParsingResult<T> : IParsingResult<T>
    {
        public SuccessfulParsingResult(T result)
        {
            Result = result;
        }

        public bool Success => true;
        public string ErrorMessage => "";
        public T Result { get; }
    }
}