using System;
using System.Linq;
using CanvasCommon.CommandLineParsing.CoreOptionTypes;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Illumina.Common.FileSystem;

namespace EvaluateCNV
{
    internal static class Program
    {
        public static int Main(string[] args)
        {
            try
            {
                var result = MainHelper(args);
                return result;
            }
            catch (Exception e)
            {
                Console.Error.WriteLine(e);
                return -1;
            }
        }

        private static int MainHelper(string[] args)
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
        private static readonly ValueOption<string> KmerFa = ValueOption<string>.Create("Kmer.fa file", "k", "kmerFa");
        private static readonly FlagOption SkipDiploid = new FlagOption("Skip diploid calls", "d", "skipDiploid");
        private static readonly FlagOption Help = new FlagOption("show this message and exit", "h", "help");

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
                KmerFa,
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
            string kmerFa = parseInput.Get(KmerFa);
            var help = parseInput.Get(Help);
            return ParsingResult<EvaluateCnvOptions>.SuccessfulResult(new EvaluateCnvOptions(baseFileName, roiBed, heterogeneityFraction,
                dqscoreThreshold, splitBySize, skipDiploid, kmerFa, help));

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
        public string KmerFa { get; }

        public bool Help { get; }

        public EvaluateCnvOptions(string baseFileName, IFileLocation roiBed, double heterogeneityFraction, double? dqscoreThreshold,
            bool splitBySize, bool skipDiploid, string kmerFa, bool help)
        {
            BaseFileName = baseFileName;
            RoiBed = roiBed;
            HeterogeneityFraction = heterogeneityFraction;
            DQscoreThreshold = dqscoreThreshold;
            SplitBySize = splitBySize;
            SkipDiploid = skipDiploid;
            KmerFa = kmerFa;
            Help = help;
        }
    }
}
