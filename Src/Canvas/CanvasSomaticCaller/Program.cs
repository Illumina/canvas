 using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using NDesk.Options;
using SequencingFiles;

namespace CanvasSomaticCaller
{
    class Program
    {
        static void ShowHelp(OptionSet p)
        {
            Console.WriteLine("Usage: CanvasSomaticCaller.exe [OPTIONS]+");
            Console.WriteLine("Make discrete-valued copy number calls for a somatic sample.");
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
            string referenceFolder = null;
            string name = "SAMPLE";
            string truthDataPath = null;
            string somaticVCFPath = null;
            bool needHelp = false;
            string bedPath = null;
            string ploidyBedPath = null;
            string ffpeOutliersPath = null;
            bool isEnrichment = false;
            bool isDbsnpVcf = false;
            double? localSDmetric = null;
            // Parameters, for parameter-sweep:
            float deviationFactor = SomaticCaller.DefaultDeviationFactor;
            double coverageWeighting = SomaticCaller.DefaultCoverageWeighting;
            int minimumCallSize = SomaticCaller.DefaultMinimumCallSize;
            double precisionWeightingFactor = SomaticCaller.DefaultPrecisionWeightingFactor;
            float? userPurity = null;
            float? userPloidy = null;

            OptionSet p = new OptionSet()
            {
                { "i|infile=", "file containing bins, their counts, and assigned segments (obtained from CanvasPartition.exe)",  v => inFile = v },
                { "v|varfile=", "file containing variant frequencies (obtained from CanvasSNV.exe)", v => variantFrequencyFile = v },
                { "o|outfile=", "file name prefix to ouput copy number calls to outfile.vcf", v => outFile = v },
                { "r|reference=", "folder that contains both genome.fa and GenomeSize.xml", v => referenceFolder = v },
                { "n|name=", "sample name for output VCF header (optional)", v => name = v },
                { "t|truth=", "path to vcf/bed with CNV truth data (optional)", v => truthDataPath = v },
                { "h|help", "show this message and exit", v => needHelp = v != null },
                { "e|enrichment", "flag indicating this is enrichment data", v => isEnrichment = v != null },
                { "s|somaticvcf=", "somatic vcf file - optionally used for purity estimation", v => somaticVCFPath = v },
                { "b|bedfile=", "bed file containing regions to exclude from calling", v => bedPath = v},
                { "p|ploidyBedFile=", "bed file specifying reference ploidy (e.g. for sex chromosomes) (optional)", v => ploidyBedPath = v},
                { "f|localSDFile=", "text file with localSD metric (calculate within CanvasClean) (optional)", v => ffpeOutliersPath = v},
                { "d|dbsnpvcf", "flag indicating a dbSNP VCF file is used to generate the variant frequency file", v => isDbsnpVcf = v != null },
                { "D|deviation=", "INTERNAL: best deviation parameter", v => deviationFactor = float.Parse(v) }, 
                { "C|coverageweight=", "INTERNAL: coverage weighting", v => coverageWeighting = float.Parse(v) }, 
                { "M|minimumcall=", "INTERNAL: minimum call size", v => minimumCallSize = int.Parse(v) }, 
                { "P|precisionweight=", "INTERNAL: precision weighting factor", v => precisionWeightingFactor = double.Parse(v) }, 
                { "u|definedpurity=", "INTERNAL: user pre-defined purity", v => userPurity = float.Parse(v) }, 
                { "l|definedploidy=", "INTERNAL: user pre-defined ploidy", v => userPloidy = float.Parse(v) }, 
            };

            List<string> extraArgs = p.Parse(args);

            if (extraArgs.Count > 0)
            {
                Console.WriteLine("Error: Argument '{0}' not understood", extraArgs[0]);
                needHelp = true;
            }

            if (needHelp)
            {
                ShowHelp(p);
                return 0;
            }

            if (inFile == null || outFile == null || referenceFolder == null)
            {
                ShowHelp(p);
                return 0;
            }

            if (!File.Exists(inFile))
            {
                Console.WriteLine("Canvas error: File {0} does not exist! Exiting.", inFile);
                return 1;
            }

            if (!File.Exists(variantFrequencyFile))
            {
                Console.WriteLine("Canvas error: File {0} does not exist! Exiting.", variantFrequencyFile);
                return 1;
            }

            if (!File.Exists(Path.Combine(referenceFolder, "GenomeSize.xml")))
            {
                Console.WriteLine("Canvas error: File {0} does not exist! Exiting.", Path.Combine(referenceFolder, "GenomeSize.xml"));
                return 1;
            }

            if (userPurity == null && userPloidy != null || userPurity != null && userPloidy == null)
            {
                Console.WriteLine("Canvas error: both definedpurity and definedploidy should be specified");
                return 1;
            }
            SomaticCaller caller = new SomaticCaller();
            // Set parameters:
            caller.TruthDataPath = truthDataPath;
            caller.SomaticVCFPath = somaticVCFPath;
            caller.IsEnrichment = isEnrichment;
            caller.IsDbsnpVcf = isDbsnpVcf;
            caller.DeviationFactor = deviationFactor;
            caller.MinimumCallSize = minimumCallSize;
            caller.CoverageWeighting = coverageWeighting;
            caller.userPurity = userPurity;
            caller.userPloidy = userPloidy;
            caller.PrecisionWeightingFactor = precisionWeightingFactor;
            if (!string.IsNullOrEmpty(ploidyBedPath))
            {
                caller.LoadReferencePloidy(ploidyBedPath);
            }
            if (!string.IsNullOrEmpty(ffpeOutliersPath))
            {
                localSDmetric = CanvasCommon.CanvasIO.ReadLocalSDFromTextFile(ffpeOutliersPath);
            }
            caller.LoadBedFile(bedPath);
            return caller.CallVariants(inFile, variantFrequencyFile, outFile, referenceFolder, name, localSDmetric);
        }
    }
}
