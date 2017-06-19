using CanvasCommon.CommandLineParsing.CoreOptionTypes;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Illumina.Common.FileSystem;

namespace Canvas.CommandLineParsing
{
    public class TumorNormalEnrichmentModeParser : ModeParser<TumorNormalEnrichmentInput>
    {
        private static readonly CommonOptionsParser CommonOptionsParser = new CommonOptionsParser();
        private static readonly SingleSampleCommonOptionsParser SingleSampleCommonOptionsParser = new SingleSampleCommonOptionsParser();
        private static readonly TumorNormalOptionsParser TumorNormalOptionsParser = new TumorNormalOptionsParser();
        private static readonly FileOption NormalBam = FileOption.CreateRequired("normal sample .bam file", "normal-bam");
        private static readonly FileOption Manifest = FileOption.CreateRequired("Nextera manifest file", "manifest");

        public TumorNormalEnrichmentModeParser(string name, string description) : base(name, description)
        {
        }

        public override ParsingResult<TumorNormalEnrichmentInput> GetResult(SuccessfulResultCollection result, CommonOptions commonOptions)
        {
            SingleSampleCommonOptions singleSampleCommonOptions = result.Get(SingleSampleCommonOptionsParser);
            TumorNormalOptions tumorNormalOptions = result.Get(TumorNormalOptionsParser);
            var normalBam = result.Get(NormalBam);
            var manifest = result.Get(Manifest);
            return ParsingResult<TumorNormalEnrichmentInput>.SuccessfulResult(
                new TumorNormalEnrichmentInput(commonOptions, singleSampleCommonOptions, tumorNormalOptions, normalBam, manifest));
        }

        public override ParsingResult<IModeRunner> GetRunner(TumorNormalEnrichmentInput result)
        {
            return ParsingResult<IModeRunner>.SuccessfulResult(new TumorNormalEnrichmentRunner(result));
        }

        public override OptionCollection<IModeLauncher> GetModeOptions()
        {
            return new OptionCollection<IModeLauncher>
            {
                TumorNormalOptionsParser, NormalBam, Manifest, CommonOptionsParser, SingleSampleCommonOptionsParser
            };
        }
    }

    public class TumorNormalEnrichmentInput
    {
        public TumorNormalEnrichmentInput(CommonOptions commonOptions, SingleSampleCommonOptions singleSampleCommonOptions, TumorNormalOptions tumorNormalOptions, IFileLocation normalBam, IFileLocation manifest)
        {
            CommonOptions = commonOptions;
            SingleSampleCommonOptions = singleSampleCommonOptions;
            TumorNormalOptions = tumorNormalOptions;
            NormalBam = normalBam;
            Manifest = manifest;
        }

        public CommonOptions CommonOptions { get; }
        public SingleSampleCommonOptions SingleSampleCommonOptions { get; }
        public TumorNormalOptions TumorNormalOptions { get; }
        public IFileLocation NormalBam { get; }
        public IFileLocation Manifest { get; }
    }
}