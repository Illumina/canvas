using System.Collections.Generic;
using Isas.SequencingFiles;

namespace CanvasCommon
{
    public class MultiSampleGenomicBin
    {
        public GenomicBin Bin { get; }
        public List<float> Counts { get; }

        public MultiSampleGenomicBin(GenomicBin genomicBin, List<float> counts)
        {
            Bin = genomicBin;
            Counts = counts;
        }
    }

    public class GenomicBin
    {
        public GenomicBin()
        {
        }

        public GenomicBin(string chromosome, Interval interval)
        {
            Chromosome = chromosome;
            Interval = interval;
        }

        public GenomicBin(string chromosome, Interval interval, int gc)
        {
            Chromosome = chromosome;
            Interval = interval;
            GC = gc;
        }

        public string Chromosome { get; set; }
        public int GC { get; set; }
        public Interval Interval { get; set; }
    }

    /// <summary>
    /// Container that holds information about a genomic interval and how many reads were observed therein.
    /// </summary>
    public class SampleGenomicBin
    {
        private readonly GenomicBin _genomicBin;

        public float Count { get; set; }
        public double CountDeviation { get; set; }

        //public SampleGenomicBin()
        //{
        //    _genomicBin = new GenomicBin();
        //    GenomicBin.Chromosome = null;
        //    GenomicBin.Interval = new Interval();
        //    GenomicBin.GC = -1;
        //    this.CountDeviation = -1;
        //}

        public SampleGenomicBin(string chr, int start, int stop, int gc, float count, double MadOfDIffs)
        {
            _genomicBin = new GenomicBin();
            GenomicBin.Chromosome = chr;
            GenomicBin.Interval = new Interval(start, stop);
            GenomicBin.GC = gc;
            this.Count = count;
            this.CountDeviation = MadOfDIffs;
        }

        public SampleGenomicBin(string chr, int start, int stop, int gc)
        {
            _genomicBin = new GenomicBin();
            GenomicBin.Chromosome = chr;
            GenomicBin.Interval = new Interval(start, stop);
            GenomicBin.GC = gc;
            this.CountDeviation = -1;
        }

        public SampleGenomicBin(string chr, int start, int stop, float count)
        {
            _genomicBin = new GenomicBin();
            GenomicBin.Chromosome = chr;
            GenomicBin.Interval = new Interval(start, stop);
            this.Count = count;
        }

        public SampleGenomicBin(string chr, int start, int stop, int gc, float count)
        {
            _genomicBin = new GenomicBin(chr, new Interval(start, stop), gc);
            this.Count = count;
            this.CountDeviation = -1;
        }

        public int Size
        {
            get { return GenomicBin.Interval.OneBasedEnd - GenomicBin.Interval.OneBasedStart; }
        }
        public int Start
        {
            get { return GenomicBin.Interval.OneBasedStart; }
            //set { GenomicBin.Interval.OneBasedStart = value; }
        }
        public int Stop
        {
            get { return GenomicBin.Interval.OneBasedEnd; }
            //set { GenomicBin.Interval.OneBasedEnd = value; }
        }

        public GenomicBin GenomicBin
        {
            get { return _genomicBin; }
        }


        public bool IsSameBin(SampleGenomicBin bin)
        {
            return GenomicBin.Chromosome == bin.GenomicBin.Chromosome && GenomicBin.Interval.OneBasedStart == bin.GenomicBin.Interval.OneBasedStart && GenomicBin.Interval.OneBasedEnd == bin.GenomicBin.Interval.OneBasedEnd;
        }
    }
}
