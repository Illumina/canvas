using CanvasCommon.CommandLineParsing.CoreOptionTypes;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Illumina.Common.FileSystem;

namespace Canvas.CommandLineParsing
{
    public class TumorNormalWgsModeParser : ModeParser<TumorNormalWgsInput>
    {
        private static readonly TumorNormalOptionsParser TumorNormalOptionsParser = new TumorNormalOptionsParser();
        private static readonly SingleSampleCommonOptionsParser SingleSampleCommonOptionsParser = new SingleSampleCommonOptionsParser();


        public TumorNormalWgsModeParser(string name, string description) : base(name, description)
        {
        }

        public override ParsingResult<TumorNormalWgsInput> GetResult(SuccessfulResultCollection result, CommonOptions commonOptions)
        {
            var singleSampleCommonOptions = result.Get(SingleSampleCommonOptionsParser);
            TumorNormalOptions tumorNormalOptions = result.Get(TumorNormalOptionsParser);
            return ParsingResult<TumorNormalWgsInput>.SuccessfulResult(
                new TumorNormalWgsInput(commonOptions, singleSampleCommonOptions, tumorNormalOptions));

        }

        public override ParsingResult<IModeRunner> GetRunner(TumorNormalWgsInput result)
        {
            return ParsingResult<IModeRunner>.SuccessfulResult(new TumorNormalWgsRunner(result));
        }

        public override OptionCollection<IModeLauncher> GetModeOptions()
        {
            return new OptionCollection<IModeLauncher>
            {
                TumorNormalOptionsParser, SingleSampleCommonOptionsParser
            };
        }
    }

    public class TumorNormalWgsInput
    {
        public TumorNormalWgsInput(CommonOptions commonOptions, SingleSampleCommonOptions singleSampleCommonOptions, TumorNormalOptions tumorNormalOptions)
        {
            CommonOptions = commonOptions;
            SingleSampleCommonOptions = singleSampleCommonOptions;
            TumorNormalOptions = tumorNormalOptions;
        }

        public CommonOptions CommonOptions { get; }
        public SingleSampleCommonOptions SingleSampleCommonOptions { get; }
        public TumorNormalOptions TumorNormalOptions { get; }
    }

    internal class TumorNormalOptionsParser : Option<TumorNormalOptions>
    {
        private static readonly FileOption TumorBam = FileOption.CreateRequired("tumor sample .bam file", "b", "bam");
        private static readonly FileOption SomaticVcf = FileOption.Create(".vcf file of somatic small variants found in the tumor sample. " +
                                                                          "Used as a fallback for estimating the reported tumor sample purity. " +
                                                                          "This option has no effect on called CNVs", "somatic-vcf");

        public override OptionCollection<TumorNormalOptions> GetOptions()
        {
            return new OptionCollection<TumorNormalOptions>
            {
               TumorBam, SomaticVcf
            };
        }

        public override ParsingResult<TumorNormalOptions> Parse(SuccessfulResultCollection parseInput)
        {
            var tumorBam = parseInput.Get(TumorBam);
            var somaticVcf = parseInput.Get(SomaticVcf);

            return ParsingResult<TumorNormalOptions>.SuccessfulResult(
                new TumorNormalOptions(tumorBam, somaticVcf));
        }
    }

    public class TumorNormalOptions
    {
        public IFileLocation TumorBam { get; }
        public IFileLocation SomaticVcf { get; }

        public TumorNormalOptions(IFileLocation tumorBam, IFileLocation somaticVcf)
        {
            TumorBam = tumorBam;
            SomaticVcf = somaticVcf;
        }
    }
}