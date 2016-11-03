using System.Collections.Generic;
using System.Linq;
using CanvasCommon.CommandLineParsing.CoreOptionTypes;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Isas.Shared.DataTypes;
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
        private static readonly FileOption PloidyVcf = FileOption.Create("multisample .vcf file containing regions of known ploidy. Copy number calls matching the known ploidy in these regions will be considered non-variant", "ploidy-bed");
        private static readonly FileOption PopulationBAlleleSites = FileOption.CreateRequired("vcf containing SNV b-allele sites in the population (only sites with PASS in the filter column will be used)", "population-b-allele-vcf");
        private static readonly FileOption SampleBAlleleSites = FileOption.CreateRequired("vcf containing SNV b-allele sites (only sites with PASS in the filter column will be used)", "b-allele-vcf");
        private static readonly FileOption CommonCnvsBed = FileOption.Create(".bed file containing regions of known common CNVs", "common-cnvs-bed");
        private static readonly MultiValueOption<string> Proband = new MultiValueOption<string>(StringOption.Create("Proband sample name", "proband"));
        private static readonly StringOption Mother = StringOption.Create("Mother sample name", "mother");
        private static readonly StringOption Father = StringOption.Create("Father sample name", "father");
        private static readonly ExclusiveFileOption BAlleleSites = ExclusiveFileOption.CreateRequired(SampleBAlleleSites, PopulationBAlleleSites);


        public override OptionCollection<SmallPedigreeOptions> GetOptions()
        {
            return new OptionCollection<SmallPedigreeOptions>
            {
                Bams, PloidyVcf, SampleNames, BAlleleSites, CommonCnvsBed, Proband, Mother, Father
            };
        }

        public override ParsingResult<SmallPedigreeOptions> Parse(SuccessfulResultCollection parseInput)
        {
            var bams = parseInput.Get(Bams);
            var sampleNames = parseInput.Get(SampleNames);
            var ploidyVcfs = parseInput.Get(PloidyVcf);
            var bAlleleSites = parseInput.Get(BAlleleSites);
            var commonCnvsBed = parseInput.Get(CommonCnvsBed);
            var pedigreeInfo = parseInput.Get(PedigreeInfo);


            List<SmallPedigreeSampleOptions> samples = new List<SmallPedigreeSampleOptions>();
            for (int sampleIndex = 0; sampleIndex < bams.Count; sampleIndex++)
                samples.Add(new SmallPedigreeSampleOptions(sampleNames[sampleIndex], bams[sampleIndex],
                    bAlleleSites[sampleIndex], ploidyVcfs[sampleIndex]));

            return ParsingResult<SmallPedigreeOptions>.SuccessfulResult(new SmallPedigreeOptions(samples, commonCnvsBed, pedigreeInfo));
        }
    }

    public class SmallPedigreeOptions
    {
        public List<SmallPedigreeSampleOptions> Samples { get; }
        public IFileLocation CommonCnvsBed { get; }
        public IFileLocation PedigreeInfo { get; }


        public SmallPedigreeOptions(List<SmallPedigreeSampleOptions> samples, IFileLocation commonCnvsBed, IFileLocation pedigreeInfo)
        {
            Samples = samples;
            CommonCnvsBed = commonCnvsBed;
            PedigreeInfo = pedigreeInfo;
        }
    }

    public class SmallPedigreeSampleOptions
    {
        public string SampleName { get; set; }
        public IFileLocation Bam { get; set; }
        public IFileLocation BAlleleSites { get; set; }
        public IFileLocation PloidyVcf { get; set; }

        public SmallPedigreeSampleOptions(string sampleName, IFileLocation bam, IFileLocation bAlleleSites, IFileLocation ploidyVcf)
        {
            SampleName = sampleName;
            Bam = bam;
            BAlleleSites = bAlleleSites;
            PloidyVcf = ploidyVcf;
        }
    }
}