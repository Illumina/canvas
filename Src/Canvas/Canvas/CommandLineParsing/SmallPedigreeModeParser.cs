using System;
using System.Collections.Generic;
using System.Linq;
using CanvasCommon.CommandLineParsing.CoreOptionTypes;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Isas.SequencingFiles;
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
        private static readonly FileOption PloidyVcf =
            FileOption.Create(
                "multisample .vcf file containing regions of known ploidy. Copy number calls matching the known ploidy in these regions will be considered non-variant",
                "ploidy-bed");
        private static readonly FileOption PopulationBAlleleSites =
            FileOption.CreateRequired(
                "vcf containing SNV b-allele sites in the population (only sites with PASS in the filter column will be used)",
                "population-b-allele-vcf");
        private static readonly FileOption SampleBAlleleSites =
            FileOption.CreateRequired(
                "vcf containing SNV b-allele sites (only sites with PASS in the filter column will be used)",
                "b-allele-vcf");
        private static readonly FileOption CommonCnvsBed =
            FileOption.Create(".bed file containing regions of known common CNVs", "common-cnvs-bed");
        private static readonly MultiValueOption<string> Proband =
            new MultiValueOption<string>(StringOption.Create("Proband sample name", "proband"));
        private static readonly StringOption Mother = StringOption.Create("Mother sample name", "mother");
        private static readonly StringOption Father = StringOption.Create("Father sample name", "father");
        private static readonly ExclusiveFileOption BAlleleSites = ExclusiveFileOption.CreateRequired(
            SampleBAlleleSites, PopulationBAlleleSites);


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
            var ploidyVcfs = parseInput.Get(PloidyVcf);
            var bAlleleSites = parseInput.Get(BAlleleSites);
            var mother = parseInput.Get(Mother);
            var father = parseInput.Get(Father);
            var proband = parseInput.Get(Proband);
            var commonCnvsBed = parseInput.Get(CommonCnvsBed);
            
            List<SmallPedigreeSampleOptions> samples = new List<SmallPedigreeSampleOptions>();
            foreach (var sample in sampleNameToBam)
                samples.Add(new SmallPedigreeSampleOptions(sample.Key, GetSampleType(sample.Key, mother, father, proband) , sample.Value));

            return
                ParsingResult<SmallPedigreeOptions>.SuccessfulResult(new SmallPedigreeOptions(samples, commonCnvsBed,
                    pedigreeInfo));
        }

        private SampleType GetSampleType(string sampleName, string mother, string father, List<string> probands)
        {
            int uniqSampleCounter = 0;
            if (sampleName == mother)
                uniqSampleCounter++;
            if (sampleName == father)
                uniqSampleCounter++;
            if (probands.Any(proband => sampleName == proband))
                uniqSampleCounter++;
            if (uniqSampleCounter > 1)
                throw new ArgumentException($"Sample {sampleName} must only have one sample type (mother | string | proband)");

            if (sampleName == mother)
                return SampleType.Mother;
            if (sampleName == father)
                return SampleType.Father;
            if (probands.Any(proband => sampleName == proband))
                return SampleType.Proband;   
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
                    if (sampleNames.Count != 1)
                        throw new ArgumentException($"Bam file '{bam}' must contain reads from one sample only");
                    map.Add(sampleNames.Single(), bam);
                }
            }
            return map;
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