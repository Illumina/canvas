using CanvasCommon.CommandLineParsing.CoreOptionTypes;
using CanvasCommon.CommandLineParsing.OptionProcessing;

namespace Canvas.CommandLineParsing
{
    public class TumorNormalEnrichmentModeParser : ModeParser
    {
        private static readonly CommonOptionsParser CommonOptionsParser = new CommonOptionsParser();
        private static readonly TumorNormalOptionsParser TumorNormalOptionsParser = new TumorNormalOptionsParser();
        private static readonly FileOption NormalBam = FileOption.CreateRequired("normal sample .bam file", "normal-bam");
        private static readonly FileOption Manifest = FileOption.CreateRequired("Nextera manifest file", "manifest");

        public TumorNormalEnrichmentModeParser(string name, string description) : base(name, description)
        {
        }

        public override ParsingResult<IModeRunner> Parse(SuccessfulResultCollection parseInput)
        {
            CommonOptions commonOptions = parseInput.Get(CommonOptionsParser);
            TumorNormalOptions tumorNormalOptions = parseInput.Get(TumorNormalOptionsParser);
            var normalBam = parseInput.Get(NormalBam);
            var manifest = parseInput.Get(Manifest);
            return ParsingResult<IModeRunner>.SuccessfulResult(new TumorNormalEnrichmentRunner(commonOptions, tumorNormalOptions, normalBam, manifest));
        }

        public override OptionCollection<IModeRunner> GetOptions()
        {
            return new OptionCollection<IModeRunner>
            {
                TumorNormalOptionsParser, NormalBam, Manifest, CommonOptionsParser
            };
        }
    }
}