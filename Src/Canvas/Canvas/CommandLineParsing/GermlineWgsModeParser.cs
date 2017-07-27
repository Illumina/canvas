using CanvasCommon.CommandLineParsing.CoreOptionTypes;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Illumina.Common.FileSystem;
using Isas.Framework.DataTypes;

namespace Canvas.CommandLineParsing
{
    public class GermlineWgsModeParser : ModeParser<GermlineWgsInput>
    {
        public static readonly FileOption Bam = FileOption.CreateRequired("sample .bam file", "b", "bam");
        private static readonly SingleSampleCommonOptionsParser SingleSampleCommonOptionsParser = new SingleSampleCommonOptionsParser();

        public GermlineWgsModeParser(string name, string description) : base(name, description)
        {
        }

        public override OptionCollection<IModeLauncher> GetModeOptions()
        {
            return new OptionCollection<IModeLauncher> { Bam, SingleSampleCommonOptionsParser };
        }

        public override IParsingResult<GermlineWgsInput> GetSerializedResult(SuccessfulResultCollection result, CommonOptions commonOptions)
        {
            var singleSampleCommonOptions = result.Get(SingleSampleCommonOptionsParser);
            var bam = result.Get(Bam);
            return ParsingResult<GermlineWgsInput>.SuccessfulResult(new GermlineWgsInput(commonOptions, singleSampleCommonOptions, bam));
        }

        public override IModeRunner GetRunner(GermlineWgsInput input)
        {
            return new GermlineWgsRunner(input);
        }
    }

    public class GermlineWgsInput
    {
        public GermlineWgsInput(CommonOptions commonOptions, SingleSampleCommonOptions singleSampleCommonOptions, IFileLocation bam)
        {
            CommonOptions = commonOptions;
            SingleSampleCommonOptions = singleSampleCommonOptions;
            Bam = bam;
        }

        public CommonOptions CommonOptions { get; }
        public SingleSampleCommonOptions SingleSampleCommonOptions { get; }
        public IFileLocation Bam { get; }
    }
}