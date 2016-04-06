using System;
using System.Diagnostics;
using Canvas.CommandLineParsing;

namespace Canvas
{
    public class Program
    {
        public static int Main(string[] args)
        {
            var version = FileVersionInfo.GetVersionInfo(typeof(MainParser).Assembly.Location).ProductVersion;
            var copyright = FileVersionInfo.GetVersionInfo(typeof(MainParser).Assembly.Location).LegalCopyright;
            var modeParser = new MainParser(version, copyright,
                new GermlineWgsModeParser("Germline-WGS", "CNV calling of a germline sample from whole genome sequencing data"),
                new SomaticEnrichmentModeParser("Somatic-Enrichment", "CNV calling of a somatic sample from targeted sequencing data"),
                new TumorNormalWgsModeParser("Somatic-WGS", "CNV calling of a somatic sample from whole genome sequencing data"),
                new TumorNormalEnrichmentModeParser("Tumor-normal-enrichment", "CNV calling of a tumor/normal pair from targeted sequencing data"));
            var result = modeParser.Parse(args, Console.Out, Console.Error);
            if (!result.Success)
                return -1;
            return result.Result.Launch();
        }
    }
}
