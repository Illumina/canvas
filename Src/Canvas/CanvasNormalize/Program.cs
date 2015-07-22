using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NDesk.Options;

namespace CanvasNormalize
{
    class Program
    {
        static int Main(string[] args)
        {
            CanvasCommon.Utilities.LogCommandLine(args);
            CanvasNormalizeParameters parameters = CanvasNormalizeParameters.ParseCommandLine(args);
            if (parameters == null) return 1;
            return CanvasNormalize.Run(parameters);
        }
    }

    internal class CanvasNormalizeParameters
    {
        #region Members
        // File containing bin counts for the tumor sample
        public string tumorBedPath = null;

        // Files containing bin counts for the normal samples
        public List<string> normalBedPaths = new List<string>();

        // Intermediate file containing weighted average bin counts over all the normal samples
        public string weightedAverageNormalBedPath = null;

        // Ouput file containing normalized bin counts
        public string outBedPath = null;

        // Manifest file
        public string manifestPath = null;

        // Bed file containing the reference ploidy
        public string ploidyBedPath = null;

        // TODO: enable this option when we have more than one normalization method
        //public CanvasCommon.CanvasCoverageMode normalizationMode = CanvasCommon.CanvasCoverageMode.TruncatedDynamicRange;
        #endregion

        public static CanvasNormalizeParameters ParseCommandLine(string[] args)
        {
            CanvasNormalizeParameters parameters = new CanvasNormalizeParameters();
            // Should I display a help message?
            bool needHelp = false;

            OptionSet p = new OptionSet()
                {
                    { "t|tumor=",         "bed file containing bin counts for the tumor sample", v => parameters.tumorBedPath = v },
                    { "n|normal=",        "bed file containing bin counts for a normal sample. Pass this option multiple times, once for each normal sample.", v => parameters.normalBedPaths.Add(v) },
                    { "o|out=",           "bed file to output containing normalized bin counts",     v => parameters.outBedPath = v },
                    { "w|weightedAverageNormal=",           "bed file to output containing normalized bin counts",     v => parameters.weightedAverageNormalBedPath = v },
                    { "f|manifest=",      "Nextera manifest file",                       v => parameters.manifestPath = v },
                    { "p|ploidyBedFile=", "bed file specifying reference ploidy (e.g. for sex chromosomes) (optional)", v => parameters.ploidyBedPath = v},
                    { "h|help",           "show this message and exit",                       v => needHelp = v != null },
                    //{ "m|mode=",          "normalization mode",                       v => parameters.normalizationMode = CanvasCommon.Utilities.ParseCanvasCoverageMode(v) },
                };

            Console.WriteLine("CanvasNormalize {0}", System.Reflection.Assembly.GetEntryAssembly().GetName().Version.ToString());

            List<string> extraArgs = p.Parse(args);

            // Check for required arguments. Display the help message if any of them are missing.
            if (string.IsNullOrEmpty(parameters.tumorBedPath))
            {
                Console.Error.WriteLine("Please specify the tumor bed file.");
                needHelp = true;
            }
            else if (!parameters.normalBedPaths.Any()) 
            {
                Console.Error.WriteLine("Please specify at least one normal bed file.");
                needHelp = true;
            }
            else if (string.IsNullOrEmpty(parameters.outBedPath))
            {
                Console.Error.WriteLine("Please specify an output file name.");
                needHelp = true;
            }

            if (needHelp)
            {
                ShowHelp(p);
                return null;
            }

            // Does the tumor bed file exist?
            if (!File.Exists(parameters.tumorBedPath))
            {
                Console.WriteLine("CanvasNormalize.exe: File {0} does not exist! Exiting.", parameters.tumorBedPath);
                return null;
            }

            // Does each of the normal bed files exist?
            foreach (string normalBedPath in parameters.normalBedPaths) 
            {
                if (!File.Exists(normalBedPath)) 
                {
                    Console.WriteLine("CanvasNormalize.exe: File {0} does not exist! Exiting.", normalBedPath);
                    return null;
                }
            }

            // Does the manifest file exist?
            if (!string.IsNullOrEmpty(parameters.manifestPath) && !File.Exists(parameters.manifestPath))
            {
                Console.WriteLine("CanvasNormalize.exe: File {0} does not exist! Exiting.", parameters.manifestPath);
                return null;
            }
            
            // Does the ploidy bed file exist?
            if (!string.IsNullOrEmpty(parameters.ploidyBedPath) && !File.Exists(parameters.ploidyBedPath)) 
            {
                Console.WriteLine("CanvasNormalize.exe: File {0} does not exist! Exiting.", parameters.ploidyBedPath);
                return null;
            }

            return parameters;
        }
        /// <summary>
        /// Command line help message.
        /// </summary>
        /// <param name="p">NDesk OptionSet</param>
        static void ShowHelp(OptionSet p)
        {
            Console.WriteLine("Usage: CanvasNormalize.exe [OPTIONS]+");
            Console.WriteLine("Normalize coverage of genomic intervals.");
            Console.WriteLine();
            Console.WriteLine("Options:");
            p.WriteOptionDescriptions(Console.Out);
        }

    }
}
