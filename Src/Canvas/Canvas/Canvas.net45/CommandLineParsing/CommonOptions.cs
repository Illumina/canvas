using System.Collections.Generic;
using Illumina.Common.FileSystem;

namespace Canvas.CommandLineParsing
{

    public class SingleSampleCommonOptions
    {
        public IFileLocation BAlleleSites { get; set; }
        public string SampleName { get; set; }
        public IFileLocation PloidyBed { get; }
        public bool IsDbSnpVcf { get; }

        /// <summary>
        /// common options for workflows with a primary sample
        /// </summary>
        /// <param name="bAlleleSites"></param>
        /// <param name="isDbSnpVcf"></param>
        /// <param name="ploidyBed"></param>
        /// <param name="sampleName"></param>
        public SingleSampleCommonOptions(IFileLocation bAlleleSites, bool isDbSnpVcf, IFileLocation ploidyBed, string sampleName)
        {
            System.Console.WriteLine("%%% SingleSampleCommonOptions: {0} {1} {2} {3}", bAlleleSites, isDbSnpVcf, ploidyBed, sampleName);
            BAlleleSites = bAlleleSites;
            IsDbSnpVcf = isDbSnpVcf;
            PloidyBed = ploidyBed;
            SampleName = sampleName;
        }
    }
    public class CommonOptions
    {
        public IDirectoryLocation OutputDirectory { get; }
        public IDirectoryLocation WholeGenomeFasta { get; }
        public IFileLocation KmerFasta { get; }
        public IFileLocation FilterBed { get; }
        public Dictionary<string, string> CustomParams { get; }
        public string StartCheckpoint { get; }
        public string StopCheckpoint { get; }


        // general common options
        /// <summary>
        /// common options for all workflows
        /// </summary>
        /// <param name="outputDirectory"></param>
        /// <param name="wholeGenomeFasta"></param>
        /// <param name="kmerFasta"></param>
        /// <param name="filterBed"></param>
        /// <param name="customParams"></param>
        /// <param name="startCheckpoint"></param>
        /// <param name="stopCheckpoint"></param>
        public CommonOptions(IDirectoryLocation outputDirectory, IDirectoryLocation wholeGenomeFasta, IFileLocation kmerFasta, IFileLocation filterBed, Dictionary<string, string> customParams, string startCheckpoint, string stopCheckpoint)
        {
            OutputDirectory = outputDirectory;
            WholeGenomeFasta = wholeGenomeFasta;
            KmerFasta = kmerFasta;
            FilterBed = filterBed;
            CustomParams = customParams;
            StartCheckpoint = startCheckpoint;
            StopCheckpoint = stopCheckpoint;
        }
    }
}