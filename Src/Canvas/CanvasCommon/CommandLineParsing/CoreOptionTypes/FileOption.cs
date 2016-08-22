using CanvasCommon.CommandLineParsing.OptionProcessing;
using Isas.Shared.Utilities.FileSystem;

namespace CanvasCommon.CommandLineParsing.CoreOptionTypes
{
    public class FileOption : ValueOption<IFileLocation>
    {
        private FileOption(ValueOptionInfo<string> info) : base(info)
        {
        }

        public static new FileOption CreateRequired(string description, params string[] names)
        {
            return new FileOption(new RequiredValueOptionInfo(description, names));
        }

        public static new FileOption Create(string description, params string[] names)
        {
            return new FileOption(new ValueOptionInfo<string>(false, description, names));
        }

        public override ParsingResult<IFileLocation> Parse(string value)
        {
            IFileLocation location = value == null ? null : new FileLocation(value);
            if (location == null || location.Exists)
                return ParsingResult<IFileLocation>.SuccessfulResult(location);
            return ParsingResult<IFileLocation>.FailedResult($"Error: {location} does not exist");
        }
    }
}