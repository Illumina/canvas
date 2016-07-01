using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NDesk.Options;

using SequencingFiles;
using Isas.Shared;
using CanvasCommon;

namespace CanvasNormalize
{
    class Program
    {
        static int Main(string[] args)
        {
            CanvasCommon.Utilities.LogCommandLine(args);
            CanvasNormalizeParameters parameters = CanvasNormalizeParameters.ParseCommandLine(args);
            if (parameters == null) return 1;

            CanvasNormalizeFactory factory = new CanvasNormalizeFactory(parameters);

            return new CanvasNormalize(factory.GetReferenceGenerator(), factory.GetRatioCalculator())
                .Run(parameters);
        }
    }

    internal class CanvasNormalizeParameters
    {
        #region Members
        // File containing bin counts for the tumor sample
        public IFileLocation tumorBedFile = null;

        // Files containing bin counts for the normal samples
        public List<IFileLocation> normalBedFiles = new List<IFileLocation>();

        // Intermediate file containing weighted average bin counts over all the normal samples
        public IFileLocation weightedAverageNormalBedFile = null;

        // Ouput file containing normalized bin counts
        public IFileLocation outBedFile = null;

        // Manifest file
        public IFileLocation manifestFile = null;

        // Bed file containing the reference ploidy
        public IFileLocation ploidyBedFile = null;

        public CanvasCommon.CanvasNormalizeMode normalizationMode = CanvasCommon.CanvasNormalizeMode.WeightedAverage;
        #endregion

        public static CanvasNormalizeParameters ParseCommandLine(string[] args)
        {
            CanvasNormalizeParameters parameters = new CanvasNormalizeParameters();
            parameters.normalizationMode = CanvasCommon.CanvasNormalizeMode.WeightedAverage;
            // Should I display a help message?
            bool needHelp = false;

            OptionSet p = new OptionSet()
                {
                    { "t|tumor=",         "bed file containing bin counts for the tumor sample", v => parameters.tumorBedFile = new FileLocation(v) },
                    { "n|normal=",        "bed file containing bin counts for a normal sample. Pass this option multiple times, once for each normal sample.", v => parameters.normalBedFiles.Add(new FileLocation(v)) },
                    { "o|out=",           "bed file to output containing normalized bin counts",     v => parameters.outBedFile = new FileLocation(v) },
                    { "w|weightedAverageNormal=",           "bed file to output containing reference bin counts",     v => parameters.weightedAverageNormalBedFile = new FileLocation(v) },
                    { "f|manifest=",      "Nextera manifest file",                       v => parameters.manifestFile = new FileLocation(v) },
                    { "p|ploidyBedFile=", "bed file specifying reference ploidy (e.g. for sex chromosomes) (optional)", v => parameters.ploidyBedFile = new FileLocation(v)},
                    { "h|help",           "show this message and exit",                       v => needHelp = v != null },
                    { "m|mode=",          "normalization mode (WeightedAverage/BestLR2/PCA). Default: " + parameters.normalizationMode, v => parameters.normalizationMode = CanvasCommon.Utilities.ParseCanvasNormalizeMode(v) },
                };

            Console.WriteLine("CanvasNormalize {0}", System.Reflection.Assembly.GetEntryAssembly().GetName().Version.ToString());

            List<string> extraArgs = p.Parse(args);

            // Check for required arguments. Display the help message if any of them are missing.
            if (parameters.tumorBedFile == null)
            {
                Console.Error.WriteLine("Please specify the tumor bed file.");
                needHelp = true;
            }
            else if (!parameters.normalBedFiles.Any()) 
            {
                Console.Error.WriteLine("Please specify at least one normal bed file.");
                needHelp = true;
            }
            else if (parameters.outBedFile == null)
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
            if (!parameters.tumorBedFile.Exists)
            {
                Console.WriteLine("CanvasNormalize.exe: File {0} does not exist! Exiting.", parameters.tumorBedFile);
                return null;
            }

            // Does each of the normal bed files exist?
            foreach (var normalBedFile in parameters.normalBedFiles) 
            {
                if (!normalBedFile.Exists) 
                {
                    Console.WriteLine("CanvasNormalize.exe: File {0} does not exist! Exiting.", normalBedFile.FullName);
                    return null;
                }
            }

            // Does the manifest file exist?
            if (parameters.manifestFile != null && !parameters.manifestFile.Exists)
            {
                Console.WriteLine("CanvasNormalize.exe: File {0} does not exist! Exiting.", parameters.manifestFile.FullName);
                return null;
            }
            
            // Does the ploidy bed file exist?
            if (parameters.ploidyBedFile != null && !parameters.ploidyBedFile.Exists) 
            {
                Console.WriteLine("CanvasNormalize.exe: File {0} does not exist! Exiting.", parameters.ploidyBedFile.FullName);
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
