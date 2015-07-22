using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EvaluateCNV
{
    class Program
    {
        static void Main(string[] args)
        {
            if (args.Length < 4)
            {
                Console.WriteLine("EvaluateCNV {0}",
                    System.Reflection.Assembly.GetExecutingAssembly().GetName().Version);
                Console.WriteLine("For more info see: http://sea-wiki.illumina.com/display/BIOINFO/EvaluateCNV");
                Console.WriteLine();
                Console.WriteLine("Usage info:");
                Console.WriteLine("EvaluateCNV $TruthSetPath $CNV.vcf $ExcludedRegionsBed $OutputPath [$RegionOfInterestBed]");
                return;
            }
            CNVChecker checker = new CNVChecker();
            string ROIBed = null;
            if (args.Length > 4) ROIBed = args[4];
            checker.Evaluate(args[0], args[1], args[2], args[3], ROIBed);
        }
    }
}
