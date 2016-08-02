using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CanvasCommon
{

    /// <summary>
    /// Container that holds information about a genomic interval and how many reads were observed therein.
    /// </summary>
    public class GenomicBin
    {

        public string Chromosome { get; set; }
        public int Start { get; set; }
        public int Stop { get; set; }
        public int GC { get; set; }
        public float Count { get; set; }
        public double MadOfDiffs{ get; set; }


        public GenomicBin()
        {
            this.Chromosome = null;
            this.Start = -1;
            this.Stop = -1;
            this.GC = -1;
            this.MadOfDiffs = -1;
            Count = -1;
        }

        public GenomicBin(string chr, int start, int stop, int gc, float count, double MadOfDIffs)
        {

            this.Chromosome = chr;
            this.Start = start;
            this.Stop = stop;
            this.GC = gc;
            this.Count = count;
            this.MadOfDiffs = MadOfDIffs;

        }

        public GenomicBin(string chr, int start, int stop, int gc, float count)
        {

            this.Chromosome = chr;
            this.Start = start;
            this.Stop = stop;
            this.GC = gc;
            this.Count = count;
            this.MadOfDiffs = -1;

        }

        public int Size
        {
            get { return Stop - Start; }
        }

        public bool IsSameBin(GenomicBin bin)
        {
            return Chromosome == bin.Chromosome && Start == bin.Start && Stop == bin.Stop;
        }
    }
}
