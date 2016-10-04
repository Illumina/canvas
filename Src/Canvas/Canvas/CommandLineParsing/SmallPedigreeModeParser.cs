using System.Collections.Generic;
using System.Linq;
using CanvasCommon.CommandLineParsing.CoreOptionTypes;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Isas.Shared.Utilities.FileSystem;

namespace Canvas.CommandLineParsing
{
    public class SmallPedigreeModeParser : ModeParser
    {
        private static readonly CommonOptionsParser CommonOptionsParser = new CommonOptionsParser();
        private static readonly SmallPedigreeOptionsParser SmallPedigreeOptionsParser = new SmallPedigreeOptionsParser();

        public SmallPedigreeModeParser(string name, string description) : base(name, description)
        {
        }

        public override ParsingResult<IModeRunner> Parse(SuccessfulResultCollection parseInput)
        {
            CommonOptions commonOptions = parseInput.Get(CommonOptionsParser);
            SmallPedigreeOptions smallPedigreeOptions = parseInput.Get(SmallPedigreeOptionsParser);
            return ParsingResult<IModeRunner>.SuccessfulResult(new SmallPedigreeWgsRunner(commonOptions, smallPedigreeOptions));
        }

        public override OptionCollection<IModeRunner> GetOptions()
        {
            return new OptionCollection<IModeRunner>
            {
                SmallPedigreeOptionsParser, CommonOptionsParser
            };
        }
    }

    internal class SmallPedigreeOptionsParser : Option<SmallPedigreeOptions>
    {
        private static readonly MultiValueOption<IFileLocation> Bams = new MultiValueOption<IFileLocation>(FileOption.Create("Bam files", "bams"));
        private static readonly MultiValueOption<string> SampleNames = new MultiValueOption<string>(ValueOption<string>.Create("Sample names", "names"));


        public override OptionCollection<SmallPedigreeOptions> GetOptions()
        {
            return new OptionCollection<SmallPedigreeOptions>
            {
                Bams, SampleNames
            };
        }

        public override ParsingResult<SmallPedigreeOptions> Parse(SuccessfulResultCollection parseInput)
        {
            var bams = parseInput.Get(Bams);
            var sampleNames = parseInput.Get(SampleNames);
            return ParsingResult<SmallPedigreeOptions>.SuccessfulResult(new SmallPedigreeOptions(bams, sampleNames));
        }
    }

    public class SmallPedigreeOptions
    {
        public IEnumerable<IFileLocation> Bams { get; }
        public IEnumerable<string> SampleNames { get; }

        public SmallPedigreeOptions(IEnumerable<IFileLocation> bams, IEnumerable<string> sampleNames)
        {
            Bams = bams;
            SampleNames = sampleNames;
        }
    }
}