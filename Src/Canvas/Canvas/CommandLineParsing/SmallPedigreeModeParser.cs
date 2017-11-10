using Canvas.SmallPedigree;
using CanvasCommon.CommandLineParsing.OptionProcessing;

namespace Canvas.CommandLineParsing
{
    internal class SmallPedigreeModeParser : ModeParser<SmallPedigreeInput>
    {
        private static readonly SmallPedigreeOptionsParser SmallPedigreeOptionsParser = new SmallPedigreeOptionsParser();

        public SmallPedigreeModeParser(string name, string description) : base(name, description)
        {
        }

        public override IParsingResult<SmallPedigreeInput> GetSerializedResult(SuccessfulResultCollection result, CommonOptions commonOptions)
        {
            var smallPedigreeOptions = result.Get(SmallPedigreeOptionsParser);
            return ParsingResult<SmallPedigreeInput>.SuccessfulResult(new SmallPedigreeInput(commonOptions, smallPedigreeOptions));
        }

        public override IModeRunner GetRunner(SmallPedigreeInput result)
        {
            return new SmallPedigreeWgsRunner(result);
        }

        public override OptionCollection<IModeLauncher> GetModeOptions()
        {
            return new OptionCollection<IModeLauncher>
            {
                SmallPedigreeOptionsParser
            };
        }
    }

    public class SmallPedigreeInput
    {
        public SmallPedigreeInput(CommonOptions commonOptions, SmallPedigreeOptions smallPedigreeOptions)
        {
            CommonOptions = commonOptions;
            SmallPedigreeOptions = smallPedigreeOptions;
        }

        public CommonOptions CommonOptions { get; }
        public SmallPedigreeOptions SmallPedigreeOptions { get; }
    }
}