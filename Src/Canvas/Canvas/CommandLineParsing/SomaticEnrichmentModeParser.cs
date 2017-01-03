using System.Collections.Generic;
using System.Linq;
using CanvasCommon.CommandLineParsing.CoreOptionTypes;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Isas.Shared.Utilities.FileSystem;

namespace Canvas.CommandLineParsing
{
    public class SomaticEnrichmentModeParser : ModeParser
    {
        private static readonly CommonOptionsParser CommonOptionsParser = new CommonOptionsParser();
        private static readonly SingleSampleCommonOptionsParser SingleSampleCommonOptionsParser = new SingleSampleCommonOptionsParser();
        private static readonly SomaticEnrichmentOptionsParser SomaticEnrichmentOptionsParser = new SomaticEnrichmentOptionsParser();

        public SomaticEnrichmentModeParser(string name, string description) : base(name, description)
        {
        }

        public override ParsingResult<IModeRunner> Parse(SuccessfulResultCollection parseInput)
        {
            CommonOptions commonOptions = parseInput.Get(CommonOptionsParser);
            SomaticEnrichmentOptions somaticEnrichmentOptions = parseInput.Get(SomaticEnrichmentOptionsParser);
            SingleSampleCommonOptions singleSampleCommonOptions = parseInput.Get(SingleSampleCommonOptionsParser);
            return ParsingResult<IModeRunner>.SuccessfulResult(new SomaticEnrichmentRunner(commonOptions, singleSampleCommonOptions, somaticEnrichmentOptions));
        }

        public override OptionCollection<IModeRunner> GetOptions()
        {
            return new OptionCollection<IModeRunner>
            {
                SomaticEnrichmentOptionsParser, CommonOptionsParser, SingleSampleCommonOptionsParser
            };
        }
    }

    internal class SomaticEnrichmentOptionsParser : Option<SomaticEnrichmentOptions>
    {
        private static readonly FileOption Bam = FileOption.CreateRequired("tumor sample .bam file", "b", "bam");
        private static readonly FileOption Manifest = FileOption.CreateRequired("Nextera manifest file", "manifest");
        private static readonly MultiValueOption<IFileLocation> ControlBams = new MultiValueOption<IFileLocation>(FileOption.Create("Bam file of an unmatched control sample", "control-bam"));
        private static readonly FileOption ControlBinned = FileOption.Create("Canvas .binned file containing precomputed control bin data to use for normalization", "control-binned");
        private static readonly ValueOption<uint?> ControlBinSize = ValueOption<uint?>.Create("bin size for control .binned file", "control-bin-size");
        private static readonly FileOption ControlPloidyBed = FileOption.Create(".bed file containing regions of known ploidy for the control .binned file", "control-ploidy-bed");

        public override OptionCollection<SomaticEnrichmentOptions> GetOptions()
        {
            return new OptionCollection<SomaticEnrichmentOptions>
            {
               Bam, Manifest, ControlBams, ControlBinned, ControlBinSize, ControlPloidyBed
            };
        }

        public override ParsingResult<SomaticEnrichmentOptions> Parse(SuccessfulResultCollection parseInput)
        {
            var bam = parseInput.Get(Bam);
            var manifest = parseInput.Get(Manifest);
            var controlBams = parseInput.Get(ControlBams);
            var controlBinned = parseInput.Get(ControlBinned);
            var controlBinSize = parseInput.Get(ControlBinSize);
            var controlPloidy = parseInput.Get(ControlPloidyBed);

            var controlBinnedBools = new List<bool> { controlBinned != null, controlBinSize != null, controlPloidy != null };
            if (controlBams.Any() && controlBinnedBools.Any(controlBinnedBool => controlBinnedBool))
                return ParsingResult<SomaticEnrichmentOptions>.FailedResult($"Error: option {ControlBams.Info.Name} cannot be combined with any of {ControlBinned.Info.Name}, {ControlBinSize.Info.Name}, {ControlPloidyBed.Info.Name} ");

            if (controlBinned != null && controlBinSize == null)
                return ParsingResult<SomaticEnrichmentOptions>.FailedResult($"Error: {ControlBinSize.Info.Name} is required when using the {ControlBinned.Info.Name} option");

            return ParsingResult<SomaticEnrichmentOptions>.SuccessfulResult(
                new SomaticEnrichmentOptions(
                    bam,
                    manifest,
                    controlBams,
                    controlBinned,
                    controlBinSize.HasValue ? (int)controlBinSize.Value : 0,
                    controlPloidy));
        }
    }

    public class SomaticEnrichmentOptions
    {
        public IFileLocation Bam { get; }
        public IFileLocation Manifest { get; }
        public IEnumerable<IFileLocation> ControlBams { get; }
        public IFileLocation ControlBinned { get; }
        public int ControlBinSize { get; }
        public IFileLocation ControlPloidy { get; }

        public SomaticEnrichmentOptions(IFileLocation bam, IFileLocation manifest, IEnumerable<IFileLocation> controlBams, IFileLocation controlBinned, int controlBinSize, IFileLocation controlPloidy)
        {
            Bam = bam;
            Manifest = manifest;
            ControlBams = controlBams;
            ControlBinned = controlBinned;
            ControlBinSize = controlBinSize;
            ControlPloidy = controlPloidy;
        }
    }
}