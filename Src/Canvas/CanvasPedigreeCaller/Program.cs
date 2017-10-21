using System;
using System.Collections.Generic;
using System.IO;
using Newtonsoft.Json;
using System.Linq;
using CanvasCommon;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Isas.Framework.Logging;
using Utilities = Isas.Framework.Utilities.Utilities;

namespace CanvasPedigreeCaller
{
    class Program
    {

        /// <summary>
        /// Command line help message.
        /// </summary>
        /// <param name="p">NDesk OptionSet</param>
        static void ShowHelp(OptionSet p)
        {
            Console.WriteLine("Usage: CanvasPedigreeCaller.exe [OPTIONS]+");
            Console.WriteLine("Make discrete-valued copy number calls in a pedigree.");
            Console.WriteLine();
            Console.WriteLine("Options:");
            p.WriteOptionDescriptions(Console.Out);
        }

        static int Main(string[] args)
        {
            CanvasCommon.Utilities.LogCommandLine(args);
            string outDir = null;
            var segmentFiles = new List<string>();
            var variantFrequencyFiles = new List<string>();
            var sampleTypesString = new List<string>();
            string ploidyBedPath = null;
            string referenceFolder = null;
            var sampleNames = new List<string>();
            bool needHelp = false;
            int? qScoreThresholdOption = null;
            int? dqScoreThresholdOption = null;
            string commonCnvsBedPath = null;
            string parameterconfigPath = Path.Combine(Utilities.GetAssemblyFolder(typeof(Program)), "PedigreeCallerParameters.json");

            var p = new OptionSet()
            {
                { "i|infile=",        "file containing bins, their counts, and assigned segments (obtained from CanvasPartition.exe)",  v => segmentFiles.Add(v) },
                { "v|varfile=",       "file containing variant frequencies (obtained from CanvasSNV.exe)",                              v => variantFrequencyFiles.Add(v) },
                { "t|sampleType=",    "sample types",                                                                                   v => sampleTypesString.Add(v) },
                { "o|outdir=",        "name of output directory",                                                                       v => outDir = v },
                { "r|reference=",     "reference genome folder that contains GenomeSize.xml",                                           v => referenceFolder = v },
                { "n|sampleName=",    "sample name for output VCF header (optional)",                                                   v => sampleNames.Add(v)},
                { "p|ploidyBed=",     "bed file specifying reference ploidy (e.g. for sex chromosomes) (optional)",                     v => ploidyBedPath = v },
                { "h|help",           "show this message and exit",                                                                     v => needHelp = v != null },
                { "q|qscore=",        $"quality filter threshold (default {CanvasPedigreeCaller.DefaultQualityFilterThreshold})",                            v => qScoreThresholdOption = int.Parse(v) },
                { "commoncnvs=",      "bed file with common CNVs (always include these intervals into segmentation results)",           v => commonCnvsBedPath = v },
                { "d|dqscore=",       $"de novo quality filter threshold (default {CanvasPedigreeCaller.DefaultDeNovoQualityFilterThreshold})",              v => dqScoreThresholdOption = int.Parse(v) },
                { "c|config=",        $"parameter configuration path (default {parameterconfigPath})",                                  v => parameterconfigPath = v}
            };

            var extraArgs = p.Parse(args);

            if (extraArgs.Count > 0)
            {
                Console.WriteLine("* Error: I don't understand the argument '{0}'", extraArgs[0]);
                needHelp = true;
            }

            if (needHelp)
            {
                ShowHelp(p);
                return 0;
            }

            if (!segmentFiles.Any() || !variantFrequencyFiles.Any() || string.IsNullOrEmpty(referenceFolder) || string.IsNullOrEmpty(outDir))
            {
                ShowHelp(p);
                return 0;
            }

            foreach (string segmentFile in segmentFiles)
            {
                if (File.Exists(segmentFile)) continue;
                Console.WriteLine($"CanvasPedigreeCaller.exe: File {segmentFile} does not exist! Exiting.");
                return 1;
            }

            foreach (string variantFrequencyFile in variantFrequencyFiles)
            {
                if (File.Exists(variantFrequencyFile)) continue;
                Console.WriteLine($"CanvasPedigreeCaller.exe: File {variantFrequencyFile} does not exist! Exiting.");
                return 1;
            }

            var sampleTypesEnum = sampleTypesString.Select(GetSampleType).ToList();

            if (!File.Exists(Path.Combine(referenceFolder, "GenomeSize.xml")))
            {
                Console.WriteLine($"CanvasPedigreeCaller.exe: File {Path.Combine(referenceFolder, "GenomeSize.xml")} does not exist! Exiting.");
                return 1;
            }

            if (!File.Exists(parameterconfigPath))
            {
                Console.WriteLine($"CanvasPedigreeCaller.exe: File {parameterconfigPath} does not exist! Exiting.");
                return 1;
            }

            if (commonCnvsBedPath != null)
            {
                if (!File.Exists(commonCnvsBedPath))
                {
                    Console.WriteLine($"CanvasPedigreeCaller.exe: File {commonCnvsBedPath} does not exist! Exiting.");
                    return 1;
                }
            }

            var parameterconfigFile = new FileLocation(parameterconfigPath);
            var callerParameters = Deserialize<PedigreeCallerParameters>(parameterconfigFile);

            int qScoreThreshold = CanvasPedigreeCaller.DefaultQualityFilterThreshold;
            if (qScoreThresholdOption != null)
            {
                qScoreThreshold = qScoreThresholdOption.Value;
                Console.WriteLine($"CanvasPedigreeCaller.exe: Using user-supplied quality score threshold {qScoreThresholdOption}.");
            }
            if (qScoreThreshold < 0 || qScoreThreshold >= callerParameters.MaxQscore)
                throw new IlluminaException($"Quality score threshold must be >= 0 and < {callerParameters.MaxQscore}");

            int dqScoreThreshold = CanvasPedigreeCaller.DefaultDeNovoQualityFilterThreshold;
            if (dqScoreThresholdOption != null)
            {

                dqScoreThreshold = dqScoreThresholdOption.Value;
                Console.WriteLine($"CanvasPedigreeCaller.exe: Using user-supplied de novo quality score threshold {qScoreThresholdOption}.");
            }
            if (dqScoreThreshold < 0 || dqScoreThreshold >= callerParameters.MaxQscore)
                throw new IlluminaException($"De novo quality score threshold must be >= 0 and < {callerParameters.MaxQscore}");

            var logger = new Logger(new[] { Console.Out }, new[] { Console.Error });
            var copyNumberLikelihoodCalculator = new CopyNumberLikelihoodCalculator(callerParameters.MaximumCopyNumber);
            IVariantCaller variantCaller = new VariantCaller(copyNumberLikelihoodCalculator, callerParameters, qScoreThreshold);
            var caller = new CanvasPedigreeCaller(logger, qScoreThreshold, dqScoreThreshold, callerParameters, copyNumberLikelihoodCalculator, variantCaller);

            return caller.CallVariants(variantFrequencyFiles, segmentFiles, outDir, ploidyBedPath, referenceFolder, sampleNames, commonCnvsBedPath, sampleTypesEnum);
        }

        private static SampleType GetSampleType(string sampleType)
        {
            if (Enum.TryParse(sampleType, out SampleType sampleTypeEnum))
                return sampleTypeEnum;
            throw new ArgumentException($"CanvasPedigreeCaller.exe: SampleType {sampleType} does not exist!");
        }

        private static T Deserialize<T>(IFileLocation path)
        {
            using (StreamReader reader = path.OpenText())
                return JsonConvert.DeserializeObject<T>(reader.ReadToEnd());
        }
    }
}
