using System;
using System.Collections.Generic;
using System.Linq;
using CanvasCommon.CommandLineParsing.CoreOptionTypes;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Illumina.Common.FileSystem;
using Isas.SequencingFiles;

namespace Canvas.CommandLineParsing
{
    internal class SmallPedigreeOptionsParser : Option<SmallPedigreeOptions>
    {
        public const string ProbandOptionName = "proband";
        public const string MotherOptionName = "mother";
        public const string FatherOptionName = "father";
        private static readonly MultiValueOption<IFileLocation> Bams = new MultiValueOption<IFileLocation>(GermlineWgsModeParser.Bam);
        private static readonly FileOption PloidyVcf = FileOption.Create("multisample .vcf file containing regions of known ploidy. Copy number calls matching the known ploidy in these regions will be considered non-variant", "ploidy-bed");
        private static readonly FileOption PopulationBAlleleSites = SingleSampleCommonOptionsParser.PopulationBAlleleSites;
        private static readonly FileOption SampleBAlleleSites = FileOption.CreateRequired("multisample .vcf file containing SNV b-allele sites (only sites with PASS in the filter column will be used)", "b-allele-vcf");
        private static readonly FileOption CommonCnvsBed = FileOption.Create(".bed file containing regions of known common CNVs", "common-cnvs-bed");
        private static readonly MultiValueOption<string> Proband = new MultiValueOption<string>(StringOption.Create("Proband sample name", ProbandOptionName));
        private static readonly StringOption Mother = StringOption.Create("Mother sample name", MotherOptionName);
        private static readonly StringOption Father = StringOption.Create("Father sample name", FatherOptionName);
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
                BamReader.WrapException(bam, reader => map.Add(reader.GetReadGroupSample(), bam));
            }
            return map;
        }
    }
}