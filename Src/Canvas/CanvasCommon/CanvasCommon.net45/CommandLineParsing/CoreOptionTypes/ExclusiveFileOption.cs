using System.Linq;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Illumina.Common.FileSystem;

namespace CanvasCommon.CommandLineParsing.CoreOptionTypes
{
    public class ExclusiveFileOptionResult
    {
        public IFileLocation Result { get; }
        public FileOption MatchedOption { get; }

        public ExclusiveFileOptionResult(IFileLocation result, FileOption matchedOption)
        {
            Result = result;
            MatchedOption = matchedOption;
        }
    }

    public class ExclusiveFileOption : Option<ExclusiveFileOptionResult>
    {
        private readonly FileOption _option1;
        private readonly FileOption _option2;
        private readonly bool _isRequired;

        public ExclusiveFileOption(FileOption option1, FileOption option2, bool isRequired)
        {
            _option1 = option1;
            _option2 = option2;
            _isRequired = isRequired;
        }

        public static ExclusiveFileOption CreateRequired(FileOption option1, FileOption option2)
        {
            return Create(option1, option2, true);
        }

        public static ExclusiveFileOption Create(FileOption option1, FileOption option2)
        {
            return Create(option1, option2, false);
        }

        private static ExclusiveFileOption Create(FileOption option1, FileOption option2, bool isRequired)
        {
            var updatedOption1Description = option1.Info.Description + $" (either this option or option {option2.Info.Name} is required)";
            var updatedOption2Description = option2.Info.Description + $" (either this option or option {option1.Info.Name} is required)";
            var updateOption1 = FileOption.Create(updatedOption1Description, option1.Info.Names.ToArray());
            var updateOption2 = FileOption.Create(updatedOption2Description, option2.Info.Names.ToArray());
            return new ExclusiveFileOption(updateOption1, updateOption2, isRequired);
        }

        public override OptionCollection<ExclusiveFileOptionResult> GetOptions()
        {
            return new OptionCollection<ExclusiveFileOptionResult>
            {
                _option1, _option2
            };
        }

        public override ParsingResult<ExclusiveFileOptionResult> Parse(SuccessfulResultCollection parseInput)
        {
            var result1 = parseInput.Get(_option1);
            var result2 = parseInput.Get(_option2);
            if (result1 != null && result2 != null)
                return ParsingResult<ExclusiveFileOptionResult>.FailedResult($"Please specify either option {_option1.Info.Name} or option {_option2.Info.Name}, but not both");
            if (_isRequired && result1 == null && result2 == null)
                return ParsingResult<ExclusiveFileOptionResult>.FailedResult($"Either option {_option1.Info.Name} or option {_option2.Info.Name} must be specified");
            IFileLocation result = null;
            FileOption option = null;
            if (result1 != null)
            {
                result = result1;
                option = _option1;
            }
            if (result2 != null)
            {
                result = result2;
                option = _option2;
            }
            return ParsingResult<ExclusiveFileOptionResult>.SuccessfulResult(new ExclusiveFileOptionResult(result, option));
        }
    }
}