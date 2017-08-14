using System;
using System.Collections.Generic;
using System.IO;
using Newtonsoft.Json;
using System.Linq;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Isas.Framework.Utilities;

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
            string ploidyBedPath = null;
            string pedigreeFile = null;
            string referenceFolder = null;
            var sampleNames = new List<string>();
            bool needHelp = false;
            int? qScoreThreshold = null;
            int? dqScoreThreshold = null;
            string commonCNVsbedPath = null;
            string parameterconfigPath = Path.Combine(Utilities.GetAssemblyFolder(typeof(Program)), "PedigreeCallerParameters.json");
            var caller = new CanvasPedigreeCaller();

            var p = new OptionSet()
            {
                { "i|infile=",        "file containing bins, their counts, and assigned segments (obtained from CanvasPartition.exe)",  v => segmentFiles.Add(v) },
                { "v|varfile=",       "file containing variant frequencies (obtained from CanvasSNV.exe)",                              v => variantFrequencyFiles.Add(v) },
                { "o|outdir=",        "name of output directory",                                                                       v => outDir = v },
                { "r|reference=",     "reference genome folder that contains GenomeSize.xml",                                           v => referenceFolder = v },
                { "n|sampleName=",    "sample name for output VCF header (optional)",                                                   v => sampleNames.Add(v)},
                { "f|pedigree=",      "relationship within pedigree (parents/proband)",                                                 v => pedigreeFile = v },
                { "p|ploidyBed=",     "bed file specifying reference ploidy (e.g. for sex chromosomes) (optional)",                     v => ploidyBedPath = v },
                { "h|help",           "show this message and exit",                                                                     v => needHelp = v != null },
                { "q|qscore=",        $"quality filter threshold (default {caller.QualityFilterThreshold})",                            v => qScoreThreshold = int.Parse(v) },
                { "commoncnvs=",      "bed file with common CNVs (always include these intervals into segmentation results)",           v => commonCNVsbedPath = v },
                { "d|dqscore=",       $"de novo quality filter threshold (default {caller.DeNovoQualityFilterThreshold})",              v => dqScoreThreshold = int.Parse(v) },
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

            if (!segmentFiles.Any()  || !variantFrequencyFiles.Any() || string.IsNullOrEmpty(referenceFolder) || string.IsNullOrEmpty(outDir))
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

            if (commonCNVsbedPath != null)
            {
                if (!File.Exists(commonCNVsbedPath))
                {
                    Console.WriteLine($"CanvasPedigreeCaller.exe: File {commonCNVsbedPath} does not exist! Exiting.");
                    return 1;
                }
            }

            var parameterconfigFile = new FileLocation(parameterconfigPath);
            caller.CallerParameters = Deserialize<PedigreeCallerParameters>(parameterconfigFile);

            if (pedigreeFile.IsNullOrEmpty())
            {
                Console.WriteLine($"CanvasPedigreeCaller.exe: pedigreeFile option is not used! Calling CNV variants without family information.");
                return caller.CallVariants(variantFrequencyFiles, segmentFiles, outDir, ploidyBedPath, referenceFolder, sampleNames, commonCNVsbedPath);
            }

            if (qScoreThreshold.HasValue & qScoreThreshold > 0 & qScoreThreshold < caller.CallerParameters.MaxQscore)
            {
                caller.QualityFilterThreshold = qScoreThreshold.Value;
                Console.WriteLine($"CanvasPedigreeCaller.exe: Using user-supplied quality score threshold {qScoreThreshold}.");
            }


            if (dqScoreThreshold.HasValue & dqScoreThreshold > 0 & dqScoreThreshold < caller.CallerParameters.MaxQscore)
            {
                caller.DeNovoQualityFilterThreshold = dqScoreThreshold.Value;
                Console.WriteLine($"CanvasPedigreeCaller.exe: Using user-supplied de novo quality score threshold {qScoreThreshold}.");
            }

            if (!File.Exists(pedigreeFile))
            {
                Console.WriteLine($"CanvasPedigreeCaller.exe: File {pedigreeFile} does not exist! Exiting.");
                return 1;
            }

            return caller.CallVariantsInPedigree(variantFrequencyFiles, segmentFiles, outDir, ploidyBedPath, referenceFolder, sampleNames, commonCNVsbedPath, pedigreeFile);
        }

        private static T Deserialize<T>(IFileLocation path)
        {
            using (StreamReader reader = path.OpenText())
                return JsonConvert.DeserializeObject<T>(reader.ReadToEnd());
        }
    }
}
