using System;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Illumina.Common;
using Illumina.Common.FileSystem;

namespace CanvasCommon.CommandLineParsing.CoreOptionTypes
{
    public class DirectoryCreatingOption : ValueOption<IDirectoryLocation>
    {
        private DirectoryCreatingOption(ValueOptionInfo<string> info) : base(info)
        {
        }

        public new static DirectoryCreatingOption CreateRequired(string description, params string[] names)
        {
            return new DirectoryCreatingOption(new RequiredValueOptionInfo(GetDescription(description), names));
        }
        public new static DirectoryCreatingOption Create(string description, params string[] names)
        {
            return new DirectoryCreatingOption(new ValueOptionInfo<string>(false, GetDescription(description), names));
        }

        private static string GetDescription(string description) => description + ". Will be created if it doesn't exist.";

        public override IParsingResult<IDirectoryLocation> Parse(string value)
        {
            IDirectoryLocation location = value == null ? null : new DirectoryLocation(value);
            if (location != null && !ActionExtensions.Try(() => location.Create(), out Exception e))
                return ParsingResult<IDirectoryLocation>.FailedResult($"Error: failed to create directory {location}. Exception: {e}");
            return ParsingResult<IDirectoryLocation>.SuccessfulResult(location);
        }
    }
}