using System;
using System.Collections.Generic;
using CanvasCommon.CommandLineParsing.CoreOptionTypes;

namespace CanvasCommon.CommandLineParsing.OptionProcessing
{
    public class SuccessfulResultCollection
    {
        private readonly IResultCollection _collection;

        public SuccessfulResultCollection(IResultCollection collection)
        {
            IParsingResult failedResult;
            if (!collection.Validate(out failedResult))
                throw new ArgumentException($"{nameof(collection)}: all parsing results must be successful. {failedResult.ErrorMessage}.");
            _collection = collection;
        }

        public SuccessfulResultCollection(IOption option, IParsingResult result) :
            this(new ResultCollection<object>(
                new Dictionary<IOption, IParsingResult>
                {
                    {
                        option, result
                    }
                }))
        {
        }

        public T Get<T>(Option<T> option)
        {
            return _collection.Get(option).Result;
        }
    }
}