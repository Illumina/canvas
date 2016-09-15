using System.Collections.Generic;
using Isas.SequencingFiles;

namespace CanvasCommon
{

    public class SingleSampleCount
    {
        public float Count { get; set; }

        public SingleSampleCount(float count)
        {
            this.Count = count;
        }
    }

    public class MultiSampleCount
    {
        public List<float> Count { get; set; }
        public MultiSampleCount(List<float> count)
        {
            this.Count = count;
        }
    }

    /// <summary>
    /// Container that holds information about a genomic interval and how many reads were observed therein.
    /// </summary>
    public class GenomicBin
    {

        public string Chromosome { get; set; }
        public int GC { get; set; }
        public GenomicInterval Interval { get; set; }
        public SingleSampleCount CountBin { get; set; }
        public MultiSampleCount CountBins { get; set; }

        public double MadOfDiffs { get; set; }


        public GenomicBin()
        {
            this.Chromosome = null;
            this.Interval = null;
            this.CountBin = null;
            this.GC = -1;
            this.MadOfDiffs = -1;
        }

        public GenomicBin(string chr, int start, int stop, int gc, float count, double MadOfDIffs)
        {
            this.Chromosome = chr;
            this.Interval = new GenomicInterval() { Start = start, End = stop };
            this.GC = gc;
            this.CountBin = new SingleSampleCount(count);
            this.MadOfDiffs = MadOfDIffs;
        }

        public GenomicBin(string chr, int start, int stop, int gc, List<float> count)
        {
            this.Chromosome = chr;
            this.Interval = new GenomicInterval() { Start = start, End = stop };
            this.GC = gc;
            this.CountBins = new MultiSampleCount(count);
            this.MadOfDiffs = -1;
        }


        public GenomicBin(string chr, int start, int stop, int gc, float count)
        {
            this.Chromosome = chr;
            this.Interval = new GenomicInterval() { Start = start, End = stop };
            this.GC = gc;
            this.CountBin = new SingleSampleCount(count);
            this.MadOfDiffs = -1;
        }

        public int Size
        {
            get { return Interval.End - Interval.Start; }
        }
        public int Start
        {
            get { return Interval.Start; }
            set { Interval.Start = value; }
        }
        public int Stop
        {
            get { return Interval.End; }
            set { Interval.End = value; }
        }


        public bool IsSameBin(GenomicBin bin)
        {
            return Chromosome == bin.Chromosome && Interval.Start == bin.Interval.Start && Interval.End == bin.Interval.End;
        }
    }
}
