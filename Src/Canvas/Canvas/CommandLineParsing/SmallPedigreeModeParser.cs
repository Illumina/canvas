using System;
using System.Collections.Generic;
using System.Linq;
using CanvasCommon.CommandLineParsing.CoreOptionTypes;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Illumina.Common.FileSystem;
using Isas.SequencingFiles;

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
        private static readonly MultiValueOption<IFileLocation> Bams = new MultiValueOption<IFileLocation>(FileOption.CreateRequired("Bam files", "bams"));
        private static readonly FileOption PloidyVcf =FileOption.Create("multisample .vcf file containing regions of known ploidy. Copy number calls matching the known ploidy in these regions will be considered non-variant","ploidy-bed");
        private static readonly FileOption PopulationBAlleleSites = SingleSampleCommonOptionsParser.PopulationBAlleleSites;
        private static readonly FileOption SampleBAlleleSites = FileOption.CreateRequired("multisample .vcf file containing SNV b-allele sites (only sites with PASS in the filter column will be used)","b-allele-vcf");
        private static readonly FileOption CommonCnvsBed = FileOption.Create(".bed file containing regions of known common CNVs", "common-cnvs-bed");
        private static readonly MultiValueOption<string> Proband = new MultiValueOption<string>(StringOption.Create("Proband sample name", "proband"));
        private static readonly StringOption Mother = StringOption.Create("Mother sample name", "mother");
        private static readonly StringOption Father = StringOption.Create("Father sample name", "father");
        private static readonly ExclusiveFileOption BAlleleSites = ExclusiveFileOption.CreateRequired(SampleBAlleleSites, PopulationBAlleleSites);

        public override OptionCollection<SmallPedigreeOptions> GetOptions()
        {
            return new OptionCollection<SmallPedigreeOptions>
            {
                Bams,
                PloidyVcf,
                BAlleleSites,
                CommonCnvsBed,
                Proband,
                Mother,
                Father
            };
        }

        public override ParsingResult<SmallPedigreeOptions> Parse(SuccessfulResultCollection parseInput)
        {
            var bams = parseInput.Get(Bams);
            var sampleNameToBam = MapSampleNameToBam(bams);
            var ploidyVcf = parseInput.Get(PloidyVcf);
            var bAlleleSites = parseInput.Get(BAlleleSites);
            var mother = parseInput.Get(Mother);
            var father = parseInput.Get(Father);
            var proband = parseInput.Get(Proband);
            var commonCnvsBed = parseInput.Get(CommonCnvsBed);

            List<SmallPedigreeSampleOptions> samples = new List<SmallPedigreeSampleOptions>();
            foreach (var sample in sampleNameToBam)
            {
                var sampleType = GetSampleType(sample.Key, mother, father, proband);
                samples.Add(new SmallPedigreeSampleOptions(sample.Key, sampleType, sample.Value));
            }

            return ParsingResult<SmallPedigreeOptions>.SuccessfulResult(new SmallPedigreeOptions(samples, commonCnvsBed, bAlleleSites.Result, bAlleleSites.MatchedOption.Equals(PopulationBAlleleSites), ploidyVcf));
        }

        private SampleType GetSampleType(string sampleName, string mother, string father, List<string> probands)
        {
            var isMother = sampleName == mother;
            var isFather = sampleName == father;
            var isProband = probands.Any(proband => sampleName == proband);
            if (new[] { isMother, isFather, isProband }.Where(item => item).Count() > 1)
                throw new ArgumentException($"Sample {sampleName} can only have one sample type (mother | father | proband)");
            if (isMother) return SampleType.Mother;
            if (isFather) return SampleType.Father;
            if (isProband) return SampleType.Proband;
            return SampleType.Other;
        }


        private Dictionary<string, IFileLocation> MapSampleNameToBam(List<IFileLocation> bams)
        {
            var map = new Dictionary<string, IFileLocation>();
            foreach (IFileLocation bam in bams)
            {
                using (var bamReader = new BamReader(bam.FullName))
                {
                    var sampleNames = bamReader.GetReadGroupSamples();
                    if (sampleNames.Count <1)
                        throw new ArgumentException($"Bam file '{bam}' must contain reads from one sample only");
                    map.Add(sampleNames.First(), bam);
                }
            }
            return map;
        }


    }

    public class SmallPedigreeOptions
    {
        public List<SmallPedigreeSampleOptions> Samples { get; }
        public IFileLocation CommonCnvsBed { get; }
        public IFileLocation BAlleleSites { get; }
        public bool IsPopulationBAlleleSites { get; }
        public IFileLocation MultiSamplePloidyVcf { get; }


        public SmallPedigreeOptions(List<SmallPedigreeSampleOptions> samples, IFileLocation commonCnvsBed, IFileLocation bAlleleSites, bool isPopulationBAlleleSites, IFileLocation multiSamplePloidyVcf)
        {
            Samples = samples;
            CommonCnvsBed = commonCnvsBed;
            BAlleleSites = bAlleleSites;
            MultiSamplePloidyVcf = multiSamplePloidyVcf;
            IsPopulationBAlleleSites = isPopulationBAlleleSites;
        }
    }

    public enum SampleType
    {
        Mother,
        Father,
        Other,
        Proband
    }

    public class SmallPedigreeSampleOptions
    {
        public string SampleName { get; }
        public SampleType SampleType { get; }
        public IFileLocation Bam { get; }

        public SmallPedigreeSampleOptions(string sampleName, SampleType sampleType, IFileLocation bam)
        {
            SampleName = sampleName;
            SampleType = sampleType;
            Bam = bam;
        }
    }
}