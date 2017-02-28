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

        public GenomicBin(string chromosome, GenomicInterval interval)
        {
            Chromosome = chromosome;
            Interval = interval;
        }

        public GenomicBin(string chromosome, GenomicInterval interval, int gc)
        {
            Chromosome = chromosome;
            Interval = interval;
            GC = gc;
        }

        public string Chromosome { get; set; }
        public int GC { get; set; }
        public GenomicInterval Interval { get; set; }
    }

    /// <summary>
    /// Container that holds information about a genomic interval and how many reads were observed therein.
    /// </summary>
    public class SampleGenomicBin
    {
        private readonly GenomicBin _genomicBin;

        public float Count { get; set; }
        public double CountDeviation { get; set; }

        public SampleGenomicBin()
        {
            _genomicBin = new GenomicBin();
            GenomicBin.Chromosome = null;
            GenomicBin.Interval = new GenomicInterval();
            GenomicBin.GC = -1;
            this.CountDeviation = -1;
        }

        public SampleGenomicBin(string chr, int start, int stop, int gc, float count, double MadOfDIffs)
        {
            _genomicBin = new GenomicBin();
            GenomicBin.Chromosome = chr;
            GenomicBin.Interval = new GenomicInterval() { Start = start, End = stop };
            GenomicBin.GC = gc;
            this.Count = count;
            this.CountDeviation = MadOfDIffs;
        }

        public SampleGenomicBin(string chr, int start, int stop, int gc)
        {
            _genomicBin = new GenomicBin();
            GenomicBin.Chromosome = chr;
            GenomicBin.Interval = new GenomicInterval() { Start = start, End = stop };
            GenomicBin.GC = gc;
            this.CountDeviation = -1;
        }


        public SampleGenomicBin(string chr, int start, int stop, int gc, float count)
        {
            _genomicBin = new GenomicBin(chr, new GenomicInterval() { Start = start, End = stop }, gc);
            this.Count = count;
            this.CountDeviation = -1;
        }

        public int Size
        {
            get { return GenomicBin.Interval.End - GenomicBin.Interval.Start; }
        }
        public int Start
        {
            get { return GenomicBin.Interval.Start; }
            set { GenomicBin.Interval.Start = value; }
        }
        public int Stop
        {
            get { return GenomicBin.Interval.End; }
            set { GenomicBin.Interval.End = value; }
        }

        public GenomicBin GenomicBin
        {
            get { return _genomicBin; }
        }


        public bool IsSameBin(SampleGenomicBin bin)
        {
            return GenomicBin.Chromosome == bin.GenomicBin.Chromosome && GenomicBin.Interval.Start == bin.GenomicBin.Interval.Start && GenomicBin.Interval.End == bin.GenomicBin.Interval.End;
        }
    }
}
