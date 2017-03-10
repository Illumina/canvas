using System;
using System.Linq;
using CanvasCommon.CommandLineParsing.CoreOptionTypes;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Illumina.Common.FileSystem;

namespace EvaluateCNV
{
    class Program
    {
        static void Main(string[] args)
        {
            EvaluateCnvOptionsParser optionsParser = new EvaluateCnvOptionsParser();
            if (args.Length < 4)
            {
                ShowHelp(optionsParser);
                return;
            }
            var parsingResult = optionsParser.Parse(args.Skip(4));

            if (!parsingResult.Success)
            {
                Console.WriteLine(parsingResult.ErrorMessage);
                ShowHelp(optionsParser);
                Environment.Exit(1);
            }
            var options = parsingResult.Result;
            if (options.Help)
            {
                ShowHelp(optionsParser);
                return;
            }

            CNVChecker checker = new CNVChecker(options.DQscoreThreshold);
            checker.Evaluate(args[0], args[1], args[2], args[3], options);
        }

        public static void ShowHelp(EvaluateCnvOptionsParser optionsParser)
        {
            Console.WriteLine("EvaluateCNV {0}",
                    System.Reflection.Assembly.GetEntryAssembly().GetName().Version);
            Console.WriteLine("For more info see: http://confluence.illumina.com/display/BIOINFO/EvaluateCNV");
            Console.WriteLine();
            Console.WriteLine("Usage info:");
            Console.WriteLine("EvaluateCNV $TruthSetPath $CNV.vcf $ExcludedRegionsBed $OutputPath [OPTIONS]+[$RegionOfInterestBed]");
            Console.WriteLine("Options:");
            optionsParser.ShowHelp(Console.Out);
        }
    }

    public class EvaluateCnvOptionsParser : Option<EvaluateCnvOptions>
    {
        private static readonly FileOption RegionOfInterestBed = FileOption.Create("Bed file containing regions of interest to report on separately", "r", "roi");
        private static readonly ValueOption<double> HeterogeneityFraction = ValueOption<double>.CreateWithDefault(1, "HeterogeneityFraction", "het");
        private static readonly ValueOption<double?> DQscoreThreshold = ValueOption<double?>.Create("DQscore threshold", "q", "dqscore");
        private static readonly FileOption PloidyBed = FileOption.Create("Bed file specifying the regions where reference ploidy is not 2", "p", "ploidy");
        private static readonly FlagOption Help = new FlagOption("show this message and exit", "h", "help");

        public override OptionCollection<EvaluateCnvOptions> GetOptions()
        {
            return new OptionCollection<EvaluateCnvOptions>
            {
                RegionOfInterestBed,
                HeterogeneityFraction,
                DQscoreThreshold,
                PloidyBed,
                Help
            };
        }

        public override ParsingResult<EvaluateCnvOptions> Parse(SuccessfulResultCollection parseInput)
        {
            IFileLocation roiBed = parseInput.Get(RegionOfInterestBed);
            double heterogeneityFraction = parseInput.Get(HeterogeneityFraction);
            double? dqscoreThreshold = parseInput.Get(DQscoreThreshold);
            IFileLocation ploidyBed = parseInput.Get(PloidyBed);
            var help = parseInput.Get(Help);
            return ParsingResult<EvaluateCnvOptions>.SuccessfulResult(new EvaluateCnvOptions(roiBed, heterogeneityFraction, dqscoreThreshold, ploidyBed, help));

        }
    }

    public class EvaluateCnvOptions
    {
        public IFileLocation RoiBed { get; }
        public double HeterogeneityFraction { get; }
        public double? DQscoreThreshold { get; }
        public IFileLocation PloidyBed { get; }
        public bool Help { get; }

        public EvaluateCnvOptions(IFileLocation roiBed, double heterogeneityFraction, double? dqscoreThreshold, IFileLocation ploidyBed, bool help)
        {
            RoiBed = roiBed;
            HeterogeneityFraction = heterogeneityFraction;
            DQscoreThreshold = dqscoreThreshold;
            PloidyBed = ploidyBed;
            Help = help;
        }
    }
}
