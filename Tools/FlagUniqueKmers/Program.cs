using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using SequencingFiles;

namespace FlagUniqueKmers
{
    class Program
    {
        static void Main(string[] args)
        {
            //CheckFlags.CheckUniqueness();
            //CheckFlags.Main(args[0], args[1]);

            KmerChecker checker = new KmerChecker();
            if (args.Length < 2)
            {
                Console.Error.WriteLine("Usage info:");
                Console.Error.WriteLine("  FlagUniqueKmers $InputFASTA $OutputFASTA");
                return;
            }
            
            checker.Main(args[0], args[1]);
        }
    }
}
