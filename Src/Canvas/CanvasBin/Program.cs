using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Text;
using System.Threading.Tasks;
using NDesk.Options;


namespace CanvasBin
{
    public class Program
    {
        public static int Main(string[] args)
        {
            CanvasCommon.Utilities.LogCommandLine(args);
            CanvasBinParameters parameters = CanvasBinParameters.ParseCommandLine(args);
            if (parameters == null) return 1;

            Console.WriteLine("{0} Parsed command-line", DateTime.Now);
            if (parameters.coverageMode == CanvasCommon.CanvasCoverageMode.Fragment)
            {
                return (new FragmentBinner(parameters)).Bin();
            }
            else
            {
                return CanvasBin.Run(parameters);
            }
        }
    }

    public class CanvasBinParameters
    {
        #region Members
        // Fasta file of the reference genome.
        // The start of unique mers in the file should be uppercase.
        // Everything else should be lowercase.
        public string referenceFile = null;

        // If specified binning will only be performed on this chromosome
        // Use for parallelization
        public string chromosome = null;

        // BAM file containing bowtie alignments.
        // The mer length should match that used when generating the fasta file.
        // Might be generated with a linux command like "cat read1.fq | cut -b 1-35 | bowtie --quiet -q -v0 -m1 -S hg19 - | samtools view -Sbo out.bam -"
        public string bamFile = null;

        // BED file to filter on.
        public string filterFile = null;

        // Desired median counts per bin.
        public int countsPerBin = -1;

        // Bin size
        public int binSize = -1;

        // Calculate bin size only
        public bool binSizeOnly = false;

        // File to write to.
        public string outFile = null;
        public List<string> intermediatePaths = new List<string>();
        public bool isPairedEnd = false;

        // coverage modes: 0 = binary (each kmer is either hit or not), 1 = dynamic range (count up to 255 hits per kmer),
        // 2 = dynamic with outlier removal (cap coverage for a kmer at 3 * median coverage across kmers in bin)
        public CanvasCommon.CanvasCoverageMode coverageMode = CanvasCommon.CanvasCoverageMode.TruncatedDynamicRange;

        // Nextera manifest file specifying the targeted regions
        public string manifestFile = null;

        public string predefinedBinsFile = null;

        #endregion

        public static CanvasBinParameters ParseCommandLine(string[] args)
        {
            CanvasBinParameters parameters = new CanvasBinParameters();
            // Should I display a help message?
            bool needHelp = false;

            OptionSet p = new OptionSet()
                {
                    { "b|bam=",           "bam file containing unique alignments", v => parameters.bamFile = v },
                    { "r|reference=",     "Canvas-ready reference fasta file", v => parameters.referenceFile = v },
                    { "c|chr=",           "for bam input, only work on this chromosome. Output intermediate binary data. Must follow-up with a single CanvasBin call passing all the intermediate binary data files (see -i option)", v => parameters.chromosome = v},
                    { "i|infile=",        "intermediate binary data file from individual chromosome. Pass this option multiple times, once for each chromosome", v => parameters.intermediatePaths.Add(v)}, 
                    { "f|filter=",        "bed file containing regions to ignore",             v => parameters.filterFile = v },
                    { "d|bindepth=",      "median counts desired in each bin",                v => parameters.countsPerBin = Convert.ToInt32(v) },
                    { "z|binsize=",       "bin size; optional",                               v => parameters.binSize = Convert.ToInt32(v) },
                    { "o|outfile=",       "text file to output containing computed bins, or if -c option was specified the intermediate binary data file to output",     v => parameters.outFile = v },
                    { "y|binsizeonly",    "calcualte bin size and exit",                      v => parameters.binSizeOnly = v != null },
                    { "h|help",           "show this message and exit",                       v => needHelp = v != null },
                    { "p|paired-end",     "input .bam is a paired-end alignment (e.g. from Isaac)", v => parameters.isPairedEnd = v != null},
                    { "m|mode=",          "coverage measurement mode",                       v => parameters.coverageMode = CanvasCommon.Utilities.ParseCanvasCoverageMode(v) },
                    { "t|manifest=",      "Nextera manifest file",                       v => parameters.manifestFile = v },
                    { "n|bins=",          "bed file containing predefined bins",              v => parameters.predefinedBinsFile = v },
                };

            Console.WriteLine("CanvasBin {0}", System.Reflection.Assembly.GetEntryAssembly().GetName().Version.ToString());

            List<string> extraArgs = p.Parse(args);

            // Check for required arguments. Display the help message if any of them are missing.
            if (string.IsNullOrEmpty(parameters.referenceFile))
            {
                Console.Error.WriteLine("Please specify the Canvas k-uniqueness reference file.");
                needHelp = true;
            }
            else if (string.IsNullOrEmpty(parameters.outFile))
            {
                Console.Error.WriteLine("Please specify an output file name.");
                needHelp = true;
            }
            else if (parameters.coverageMode != CanvasCommon.CanvasCoverageMode.Fragment && parameters.countsPerBin == -1)
            {
                Console.Error.WriteLine("Please specify counts per bin.");
                needHelp = true;
            }
            else if (parameters.coverageMode != CanvasCommon.CanvasCoverageMode.Fragment
                && string.IsNullOrEmpty(parameters.chromosome) && parameters.intermediatePaths.Count == 0)
            {
                Console.Error.WriteLine("Please specify chromsome to measure coverage for.");
                needHelp = true;
            }

            if (needHelp)
            {
                ShowHelp(p);
                return null;
            }

            // Does the reference file exist?
            if (!File.Exists(parameters.referenceFile))
            {
                Console.WriteLine("CanvasBin.exe: File {0} does not exist! Exiting.", parameters.referenceFile);
                return null;
            }

            // Does the bam file exist?
            if (parameters.bamFile == null || !File.Exists(parameters.bamFile))
            {
                Console.WriteLine("CanvasBin.exe: Alignment input does not exist! Exiting.");
                return null;
            }

            // Does the BED file exist?
            if ((parameters.filterFile != null) && (!File.Exists(parameters.filterFile)))
            {
                Console.WriteLine("CanvasBin.exe: File {0} does not exist! Exiting.", parameters.filterFile);
                return null;
            }

            // Did the user supply a non-negative number?
            if (parameters.coverageMode != CanvasCommon.CanvasCoverageMode.Fragment && parameters.countsPerBin < 1)
            {
                Console.WriteLine("CanvasBin.exe: Median counts must be strictly positive. Exiting.");
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
            Console.WriteLine("Usage: CanvasBin.exe [OPTIONS]+");
            Console.WriteLine("Bin alignments into variable-sized genomic intervals.");
            Console.WriteLine();
            Console.WriteLine("Options:");
            p.WriteOptionDescriptions(Console.Out);

        }

    }
}
