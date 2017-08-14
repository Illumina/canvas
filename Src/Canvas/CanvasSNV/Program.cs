using System;
using System.IO;
using Illumina.Common;

namespace CanvasSNV
{
    class Program
    {
        /// <summary>
        /// Command line help message.
        /// </summary>
        /// <param name="p">NDesk OptionSet</param>
        static void ShowHelp(OptionSet p)
        {
            Console.WriteLine("Usage: CanvasSNV.exe [OPTIONS]+");
            Console.WriteLine("Parses bam file to derive allele counts.");
            Console.WriteLine();
            Console.WriteLine("Options:");
            p.WriteOptionDescriptions(Console.Out);
        }

        static int Main(string[] args)
        {
            CanvasCommon.Utilities.LogCommandLine(args);
            string chromosome = null; 
            string vcfPath = null;
            string bamPath = null;
            string outputPath = null;
            string sampleName = null;
            bool isSomatic = false;
            bool isDbSnpVcf = false; // assume vcf file comes from a small variant caller (Strelka)
            int minMapQ = 0; // only use reads with MAPQ greater than this number
            bool needHelp = false;

            var p = new OptionSet()
            {
                { "c|chromosome=",     "chromosome namne",                                                                               v => chromosome = v },
                { "v|vcfPath=",        "file containing small variants",                                                                 v => vcfPath = v },
                { "b|bamPath=",        "bam file",                                                                                       v => bamPath = v },
                { "o|outputPath=",     "name of output directory",                                                                       v => outputPath = v },
                { "n|sampleName=",     "sample name for output VCF header (optional)",                                                   v => sampleName = v},
                { "i|isDbSnpVcf=",     "flag to specify if vcf file contains dbSNP variants (optional)",                                 v => isDbSnpVcf = v != null },
                { "q|minMapQ=",        "mapQ threshold for vcf file (optional)",                                                         v => minMapQ = int.Parse(v)},
                { "s|isSomatic",       "flag to specify if Canvas workflow is somatic (optional)",                                       v => isSomatic = v != null },
                { "h|help",            "show this message and exit",                                                                     v => needHelp = v != null },
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

            if (string.IsNullOrEmpty(chromosome) || string.IsNullOrEmpty(outputPath))
            {
                ShowHelp(p);
                return 0;
            }

            if (!File.Exists(vcfPath))
            {
                Console.WriteLine($"CanvasSNV.exe: File {vcfPath} does not exist! Exiting.");
                return 1;
            }

            if (!File.Exists(bamPath))
            {
                Console.WriteLine($"CanvasSNV.exe: File {bamPath} does not exist! Exiting.");
                return 1;
            }

            // Handle some special cases, if the "chromosome" is a special string:
            switch (chromosome.ToLowerInvariant())
            {
                case "histogram":
                    return BuildEmpiricalHistograms(vcfPath, bamPath, outputPath);
                case "regionhistogram":
                    return BuildRegionHistograms(vcfPath, bamPath, outputPath);
                default:
                    // Ordinary chromosome name.
                    break;
            }

            // Standard logic: Process one chromosome, write output to the specified file path:
            SNVReviewer processor = new SNVReviewer(chromosome, vcfPath, bamPath, outputPath, sampleName, isDbSnpVcf, minMapQ, isSomatic);
            return processor.Run();
        }

        static int BuildEmpiricalHistograms(string oracleVCFPath, string empiricalVariantFrequencyFolder, string outputPath)
        {
            HistogramVF builder = new HistogramVF();
            return builder.BuildHistogramByCN(oracleVCFPath, empiricalVariantFrequencyFolder, outputPath);
        }

        static int BuildRegionHistograms(string oracleVCFPath, string empiricalVariantFrequencyFolder, string outputPath)
        {
            HistogramVF builder = new HistogramVF();
            return builder.SummarizeStatsByRegion(oracleVCFPath, empiricalVariantFrequencyFolder, outputPath);
        }

    }
}
