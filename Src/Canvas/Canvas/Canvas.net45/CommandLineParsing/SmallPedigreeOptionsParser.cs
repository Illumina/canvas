using System;
using System.Collections.Generic;
using System.Linq;
using CanvasCommon.CommandLineParsing.CoreOptionTypes;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Illumina.Common.FileSystem;
using Isas.SequencingFiles;
using Illumina.Common;

namespace Canvas.CommandLineParsing
{
    internal class SmallPedigreeOptionsParser : Option<SmallPedigreeOptions>
    {
        public const string ProbandOptionName = "proband";
        public const string MotherOptionName = "mother";
        public const string FatherOptionName = "father";
        public const string PloidyVcfOptionName = "ploidy-vcf";
        internal static readonly ValueOption<SampleType> SampleType = ValueOption<SampleType>.CreateWithDefault(CommandLineParsing.SampleType.Other, "Pedigree member type (either proband, mother, father or other). Default is other", "pedigree-member");
        private static readonly StringOption SampleName = StringOption.Create("sample name. Default is SM tag in RG header of the .bam", SingleSampleCommonOptionsParser.SampleName.Info.Names.ToArray());
        internal static readonly PositionalOption<IFileLocation, SampleType, string, SmallPedigreeSampleOptions> Bams =
            new PositionalOption<IFileLocation, SampleType, string, SmallPedigreeSampleOptions>(Parse, true,
                GermlineWgsModeParser.Bam,
                SampleType,
                SampleName,
                GermlineWgsModeParser.Bam.Info.Names.ToArray());
        private static readonly FileOption PloidyVcf = FileOption.Create("multisample .vcf file containing regions of known ploidy. Copy number calls matching the known ploidy in these regions will be considered non-variant", PloidyVcfOptionName);
        private static readonly FileOption PopulationBAlleleSites = SingleSampleCommonOptionsParser.PopulationBAlleleSites;
        private static readonly FileOption SampleBAlleleSites = FileOption.CreateRequired("multisample .vcf file containing SNV b-allele sites (only sites with PASS in the filter column will be used)", SingleSampleCommonOptionsParser.SampleBAlleleVcfOptionName);
        private static readonly FileOption CommonCnvsBed = FileOption.Create(".bed file containing regions of known common CNVs", "common-cnvs-bed");
        private static readonly ExclusiveFileOption BAlleleSites = ExclusiveFileOption.CreateRequired(SampleBAlleleSites, PopulationBAlleleSites);
        
        private static ParsingResult<SmallPedigreeSampleOptions> Parse(IFileLocation bam, SampleType sampleType, string sampleName)
        {
            if (sampleName == null)
            {
                Action a = () =>
                {
                    BamReader.WrapException(bam, reader =>
                    {
                        sampleName = reader.GetReadGroupSample();
                    });
                };
                if (!a.Try(out Exception e))
                    return ParsingResult<SmallPedigreeSampleOptions>.FailedResult(e.Message);
            }
            return ParsingResult<SmallPedigreeSampleOptions>.SuccessfulResult(new SmallPedigreeSampleOptions(sampleName, sampleType, bam));
        }

        public override OptionCollection<SmallPedigreeOptions> GetOptions()
        {
            return new OptionCollection<SmallPedigreeOptions>
            {
                Bams,
                PloidyVcf,
                BAlleleSites,
                CommonCnvsBed
            };
        }

        public override ParsingResult<SmallPedigreeOptions> Parse(SuccessfulResultCollection parseInput)
        {
            var bams = parseInput.Get(Bams);
            var ploidyVcf = parseInput.Get(PloidyVcf);
            var bAlleleSites = parseInput.Get(BAlleleSites);
            var commonCnvsBed = parseInput.Get(CommonCnvsBed);

            ParsingResult<SmallPedigreeOptions> failedResult;
            if (HasMoreThanOneSameSampleType(bams, out failedResult))
                return failedResult;
            return ParsingResult<SmallPedigreeOptions>.SuccessfulResult(new SmallPedigreeOptions(bams, commonCnvsBed, bAlleleSites.Result, bAlleleSites.MatchedOption.Equals(PopulationBAlleleSites), ploidyVcf));
        }

        private bool HasMoreThanOneSameSampleType(List<SmallPedigreeSampleOptions> bams, out ParsingResult<SmallPedigreeOptions> failedResult)
        {
            failedResult = null;
            if (HasMoreThanOne(bams, CommandLineParsing.SampleType.Mother, out failedResult) ||
                HasMoreThanOne(bams, CommandLineParsing.SampleType.Father, out failedResult) ||
                HasMoreThanOne(bams, CommandLineParsing.SampleType.Proband, out failedResult))
                return false;
            return true;
        }

        private bool HasMoreThanOne(List<SmallPedigreeSampleOptions> bams, SampleType sampleType, out ParsingResult<SmallPedigreeOptions> failedResult)
        {
            failedResult = null;
            var sameType = bams.Where(bam => bam.SampleType == sampleType);
            if (sameType.Count() <= 1) return false;

            var bamsSameType = string.Join(",", sameType.Select(bam => $"'{bam.Bam}'"));
            failedResult = ParsingResult<SmallPedigreeOptions>.FailedResult($"Pedigree can have at most one sample of type '{sampleType}'. Samples with these bams have the same type: {bamsSameType}");
            return true;
        }
    }
}