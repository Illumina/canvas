using System;
using System.Collections.Generic;
using System.IO;
using Isas.Shared.Utilities;
using Isas.Shared.Utilities.FileSystem;
using NDesk.Options;
using Newtonsoft.Json;

namespace CanvasDiploidCaller
{
    class Program
    {
        /// <summary>
        /// Command line help message.
        /// </summary>
        /// <param name="p">NDesk OptionSet</param>
        static void ShowHelp(OptionSet p)
        {
            Console.WriteLine("Usage: CanvasDiploidCaller.exe [OPTIONS]+");
            Console.WriteLine("Make discrete-valued copy number calls assuming a diploid baseline.");
            Console.WriteLine();
            Console.WriteLine("Options:");
            p.WriteOptionDescriptions(Console.Out);
        }

        static int Main(string[] args)
        {
            CanvasCommon.Utilities.LogCommandLine(args);
            string inFile = null;
            string outFile = null;
            string variantFrequencyFile = null;
            string ploidyBedPath = null;
            string referenceFolder = null;
            string sampleName = "SAMPLE";
            bool isDbsnpVcf = false;
            bool needHelp = false;
            string truthDataPath = null;

            string qualityScoreConfigPath = Path.Combine(Utilities.GetAssemblyFolder(typeof(Program)), "QualityScoreParameters.json");

            var p = new OptionSet()
            {
                { "i|infile=",        "file containing bins, their counts, and assigned segments (obtained from CanvasPartition.exe)",  v => inFile = v },
                { "v|varfile=",       "file containing variant frequencies (obtained from CanvasSNV.exe)",                              v => variantFrequencyFile = v },
                { "o|outfile=",       "file name prefix to ouput copy number calls to outfile.vcf",                                     v => outFile = v },
                { "r|reference=",     "reference genome folder that contains GenomeSize.xml",                                           v => referenceFolder = v },
                { "n|sampleName=",    "sample name for output VCF header (optional)",                                                   v => sampleName = v },
                { "p|ploidyBed=",     "bed file specifying reference ploidy (e.g. for sex chromosomes) (optional)",                     v => ploidyBedPath = v },
                { "d|dbsnpvcf", "flag indicating a dbSNP VCF file is used to generate the variant frequency file",                      v => isDbsnpVcf = v != null },
                { "h|help",           "show this message and exit",                                                                     v => needHelp = v != null },
                { "s|qscoreconfig=", $"parameter configuration path (default {qualityScoreConfigPath})", v => qualityScoreConfigPath = v },
                { "t|truth=", "path to vcf/bed with CNV truth data (optional)", v => truthDataPath = v },
            };

            List<string> extraArgs = p.Parse(args);
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

            if (inFile == null || outFile == null || string.IsNullOrEmpty(variantFrequencyFile) || string.IsNullOrEmpty(referenceFolder))
            {
                ShowHelp(p);
                return 0;
            }

            if (!File.Exists(inFile))
            {
                Console.WriteLine("CanvasDiploidCaller.exe: File {0} does not exist! Exiting.", inFile);
                return 1;
            }

            if (!File.Exists(variantFrequencyFile))
            {
                Console.WriteLine("Canvas error: File {0} does not exist! Exiting.", variantFrequencyFile);
                return 1;
            }

            if (!File.Exists(Path.Combine(referenceFolder, "GenomeSize.xml")))
            {
                Console.WriteLine("CanvasDiploidCaller.exe: File {0} does not exist! Exiting.", Path.Combine(referenceFolder, "GenomeSize.xml"));
                return 1;
            }

            FileLocation qscoreConfigFile = new FileLocation(qualityScoreConfigPath);
            CanvasCommon.QualityScoreParameters qscoreParametersJSON = Deserialize<CanvasCommon.QualityScoreParameters>(qscoreConfigFile);

            // Set parameters:
            CanvasDiploidCaller caller = new CanvasDiploidCaller();
            caller.IsDbsnpVcf = isDbsnpVcf;
            caller.germlineScoreParameters = qscoreParametersJSON;
            return caller.CallVariants(variantFrequencyFile, inFile, outFile, ploidyBedPath, referenceFolder, sampleName, truthDataPath);
        }
        private static T Deserialize<T>(IFileLocation path)
        {
            using (StreamReader reader = path.OpenText())
                return JsonConvert.DeserializeObject<T>(reader.ReadToEnd());
        }
    }
}
