using System;
using Canvas.CommandLineParsing;
using CanvasCommon;

namespace Canvas
{
    public class Program
    {
        public static int Main(string[] args)
        {
            try
            {
                var mainParser = new MainParser(CanvasVersionInfo.VersionString, CanvasVersionInfo.CopyrightString,
                    new GermlineWgsModeParser("Germline-WGS",
                        "CNV calling of a germline sample from whole genome sequencing data"),
                    new SomaticEnrichmentModeParser("Somatic-Enrichment",
                        "CNV calling of a somatic sample from targeted sequencing data"),
                    new TumorNormalWgsModeParser("Somatic-WGS",
                        "CNV calling of a somatic sample from whole genome sequencing data"),
                    new TumorNormalEnrichmentModeParser("Tumor-normal-enrichment",
                        "CNV calling of a tumor/normal pair from targeted sequencing data"),
                    new SmallPedigreeModeParser("SmallPedigree-WGS",
                        "CNV calling of a small pedigree from whole genome sequencing data"));
                var result = mainParser.Run(args, Console.Out, Console.Error);
                return result;
            }
            catch (Exception e)
            {
                Console.Error.WriteLine(e);
                return -1;
            }
        }
    }
}
