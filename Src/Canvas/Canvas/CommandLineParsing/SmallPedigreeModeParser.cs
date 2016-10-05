using System.Collections.Generic;
using System.Linq;
using CanvasCommon.CommandLineParsing.CoreOptionTypes;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Isas.Shared.Utilities.FileSystem;

namespace Canvas.CommandLineParsing
{
    public class SmallPedigreeModeParser : ModeParser
    {
        private static readonly CommonMultiSampleOptionsParser CommonOptionsParser = new CommonMultiSampleOptionsParser();
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
        private static readonly MultiValueOption<IFileLocation> Bams = new MultiValueOption<IFileLocation>(FileOption.CreateRequired("Bam files", "bams"));
        private static readonly MultiValueOption<IFileLocation> PloidyBed = new MultiValueOption<IFileLocation>(FileOption.CreateRequired(".bed file containing regions of known ploidy in the sample. Copy number calls matching the known ploidy in these regions will be considered non-variant", "ploidy-bed"));
        private static readonly MultiValueOption<string> SampleNames = new MultiValueOption<string>(StringOption.CreateRequired("Sample names", "names"));
        private static readonly MultiValueOption<IFileLocation> BAlleleSites = new MultiValueOption<IFileLocation>(FileOption.CreateRequired("vcf containing SNV b-allele sites (only sites with PASS in the filter column will be used)", "b-allele-vcf"));


        public override OptionCollection<SmallPedigreeOptions> GetOptions()
        {
            return new OptionCollection<SmallPedigreeOptions>
            {
                Bams, PloidyBed, SampleNames, BAlleleSites
            };
        }

        public override ParsingResult<SmallPedigreeOptions> Parse(SuccessfulResultCollection parseInput)
        {
            var bams = parseInput.Get(Bams);
            var sampleNames = parseInput.Get(SampleNames);
            var ploidyBed = parseInput.Get(PloidyBed);
            var bAlleleSites = parseInput.Get(BAlleleSites);
            return ParsingResult<SmallPedigreeOptions>.SuccessfulResult(new SmallPedigreeOptions(bams, sampleNames, ploidyBed, bAlleleSites));
        }
    }

    public class SmallPedigreeOptions
    {
        public IEnumerable<IFileLocation> Bams { get; }
        public IEnumerable<IFileLocation> BAlleleSites { get; }
        public IEnumerable<IFileLocation> PloidyBed { get; }

        public IEnumerable<string> SampleNames { get; }

        public SmallPedigreeOptions(IEnumerable<IFileLocation> bams, IEnumerable<string> sampleNames, IEnumerable<IFileLocation> ploidyBed, IEnumerable<IFileLocation> bAlleleSites)
        {
            Bams = bams;
            SampleNames = sampleNames;
            BAlleleSites = bAlleleSites;
            PloidyBed = ploidyBed;
        }
    }
}