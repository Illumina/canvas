using CanvasCommon.CommandLineParsing.CoreOptionTypes;
using CanvasCommon.CommandLineParsing.OptionProcessing;

namespace Canvas.CommandLineParsing
{
    public class CommonOptionsParser : Option<CommonOptions>
    {
        private static readonly FileOption KmerFasta = FileOption.CreateRequired("Canvas-ready reference fasta file", "r", "reference");
        private static readonly DirectoryOption Output = DirectoryOption.CreateRequired("output directory", "o", "output");
        private static readonly DirectoryOption WholeGenomeFasta = DirectoryOption.CreateRequired("folder that contains both genome.fa and GenomeSize.xml", "g", "genome-folder");
        private static readonly FileOption FilterFile = FileOption.CreateRequired(".bed file of regions to skip", "f", "filter-bed");
        private static readonly DictionaryOption CustomParameters = DictionaryOption.Create("use custom parameter for command-line tool. VALUE must contain the tool name followed by a comma and then the custom parameters.", "custom-parameters");
        private static readonly StringOption StartCheckpoint = StringOption.Create("continue analysis starting at the specified checkpoint. VALUE can be the checkpoint name or number", "c");
        private static readonly StringOption StopCheckpoint = StringOption.Create("stop analysis after the specified checkpoint is complete. VALUE can be the checkpoint name or number", "s");

        public override OptionCollection<CommonOptions> GetOptions()
        {
            return new OptionCollection<CommonOptions>
            {
                Output, KmerFasta, WholeGenomeFasta, FilterFile, CustomParameters, StartCheckpoint, StopCheckpoint
            };
        }

        public override ParsingResult<CommonOptions> Parse(SuccessfulResultCollection parseInput)
        {
            var output = parseInput.Get(Output);
            var wholeGenomeFasta = parseInput.Get(WholeGenomeFasta);
            var kmerFasta = parseInput.Get(KmerFasta);
            var filterBed = parseInput.Get(FilterFile);
            var customParameters = parseInput.Get(CustomParameters);
            string startCheckpoint = parseInput.Get(StartCheckpoint);
            string stopCheckpoint = parseInput.Get(StopCheckpoint);
            return ParsingResult<CommonOptions>.SuccessfulResult(
                new CommonOptions(
                    output,
                    wholeGenomeFasta,
                    kmerFasta,
                    filterBed,
                    customParameters,
                    startCheckpoint,
                    stopCheckpoint));
        }
    }
}