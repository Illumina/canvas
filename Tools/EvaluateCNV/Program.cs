using System;
using System.Collections.Generic;
using System.Linq;
using CanvasCommon;
using CanvasCommon.CommandLineParsing.CoreOptionTypes;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Illumina.Common.FileSystem;
using Illumina.SecondaryAnalysis.VariantCalling;

namespace EvaluateCNV
{
    internal static class Program
    {
        public static int Main(string[] args)
        {
            EvaluateCnvOptionsParser optionsParser = new EvaluateCnvOptionsParser();
            if (args.Length < 4)
            {
                ShowHelp(optionsParser, Console.Error);
                return 1;
            }
            var parsingResult = optionsParser.Parse(args.Skip(4));

            if (!parsingResult.Success)
            {
                Console.Error.WriteLine(parsingResult.ErrorMessage);
                ShowHelp(optionsParser, Console.Error);
                return 1;
            }
            var options = parsingResult.Result;
            if (options.Help)
            {
                ShowHelp(optionsParser, Console.Out);
                return 0;
            }

            CNVChecker.Evaluate(args[0], args[1], args[2], args[3], options);
            return 0;
        }

        public static void ShowHelp(EvaluateCnvOptionsParser optionsParser, System.IO.TextWriter writer)
        {
            writer.WriteLine("EvaluateCNV {0}",
                    System.Reflection.Assembly.GetEntryAssembly().GetName().Version);
            writer.WriteLine("For more info see: http://confluence.illumina.com/display/BIOINFO/EvaluateCNV");
            writer.WriteLine();
            writer.WriteLine("Usage info:");
            writer.WriteLine("EvaluateCNV $TruthSetPath $CNV.vcf $ExcludedRegionsBed $OutputDir [OPTIONS]+ [$RegionOfInterestBed]");
            writer.WriteLine("Options:");
            optionsParser.ShowHelp(writer.WriteLine);
        }
    }

    public class EvaluateCnvOptionsParser : Option<EvaluateCnvOptions>
    {
        private static readonly ValueOption<string> BaseFileName = ValueOption<string>.CreateWithDefault("EvaluateCNVResults", "Base file name (without extension)", "f");
        private static readonly FileOption RegionOfInterestBed = FileOption.Create("Bed file containing regions of interest to report on separately", "r", "roi");
        private static readonly ValueOption<double> HeterogeneityFraction = ValueOption<double>.CreateWithDefault(1, "HeterogeneityFraction", "het");
        private static readonly ValueOption<double?> DQscoreThreshold = ValueOption<double?>.Create("DQscore threshold", "q", "dqscore");
        private static readonly FlagOption SplitBySize = new FlagOption("Split by variant size", "s", "splitBySize");
        private static readonly FlagOption SkipDiploid = new FlagOption("Skip diploid calls", "d", "skipDiploid");
        private static readonly ValueOption<int> MinEntrySize = ValueOption<int>.CreateWithDefault(10000, "Minimum entry size to consider from either the truth or query files. Entries in those files that span fewer bases than this will be excluded from evaluation.", "min-size");
        private static readonly ValueOption<int> PloidyX = ValueOption<int>.CreateRequired("Reference ploidy for chromosome X (integer)", "ploidy-x");
        private static readonly ValueOption<int> PloidyY = ValueOption<int>.CreateRequired("Reference ploidy for chromosome Y (integer)", "ploidy-y");
        private static readonly ValueOption<IFileLocation> ParBed = FileOption.CreateRequired("Path to PAR bed file containing X chromosome PAR regions. Assumes chromosome Y PAR regions have been N-masked in the reference", "par-bed");
        private static readonly SinglePositionalOption<int, int, IFileLocation, (SexPloidyInfo, IFileLocation)> PloidyOption = new SinglePositionalOption<int, int, IFileLocation, (SexPloidyInfo, IFileLocation)>(Parse, false, "Instead of relying on the GT field in the query vcf to determine reference ploidy, use specified ploidy for allosomes taking into account chr X PAR regions in the provided bed file", PloidyX, PloidyY, ParBed, "ploidy");
        private static readonly FlagOption Help = new FlagOption("show this message and exit", "h", "help");

        private static IParsingResult<(SexPloidyInfo, IFileLocation)> Parse(int ploidyX, int ploidyY, IFileLocation parBed)
        {
            return ParsingResult<(SexPloidyInfo, IFileLocation)>.SuccessfulResult(
                (new SexPloidyInfo(ploidyX, ploidyY), parBed));
        }

        public override OptionCollection<EvaluateCnvOptions> GetOptions()
        {
            return new OptionCollection<EvaluateCnvOptions>
            {
                BaseFileName,
                RegionOfInterestBed,
                HeterogeneityFraction,
                DQscoreThreshold,
                SplitBySize,
                SkipDiploid,
                MinEntrySize,
                PloidyOption,
                Help
            };
        }

        public override IParsingResult<EvaluateCnvOptions> Parse(SuccessfulResultCollection parseInput)
        {
            string baseFileName = parseInput.Get(BaseFileName);
            IFileLocation roiBed = parseInput.Get(RegionOfInterestBed);
            double heterogeneityFraction = parseInput.Get(HeterogeneityFraction);
            double? dqscoreThreshold = parseInput.Get(DQscoreThreshold);
            bool splitBySize = parseInput.Get(SplitBySize);
            bool skipDiploid = parseInput.Get(SkipDiploid);
            int minEntrySize = parseInput.Get(MinEntrySize);
            var ploidyInfo = parseInput.Get(PloidyOption);
            var help = parseInput.Get(Help);
            return ParsingResult<EvaluateCnvOptions>.SuccessfulResult(new EvaluateCnvOptions(baseFileName, roiBed, heterogeneityFraction,
                dqscoreThreshold, splitBySize, skipDiploid, minEntrySize, ploidyInfo, help));

        }
    }

    public class EvaluateCnvOptions
    {
        public string BaseFileName { get; }
        public IFileLocation RoiBed { get; }
        public double HeterogeneityFraction { get; }
        public double? DQscoreThreshold { get; }
        public bool SplitBySize { get; }
        public bool SkipDiploid { get; }
        public int MinEntrySize { get; }
        public (SexPloidyInfo SexPloidyInfo, IFileLocation ParBed) PloidyInfo { get; }
        public bool Help { get; }

        public EvaluateCnvOptions(string baseFileName, IFileLocation roiBed, double heterogeneityFraction,
            double? dqscoreThreshold,
            bool splitBySize, bool skipDiploid, int minEntrySize, (SexPloidyInfo, IFileLocation) ploidyInfo,
            bool help)
        {
            BaseFileName = baseFileName;
            RoiBed = roiBed;
            HeterogeneityFraction = heterogeneityFraction;
            DQscoreThreshold = dqscoreThreshold;
            SplitBySize = splitBySize;
            SkipDiploid = skipDiploid;
            MinEntrySize = minEntrySize;
            PloidyInfo = ploidyInfo;
            Help = help;
        }
    }
}
