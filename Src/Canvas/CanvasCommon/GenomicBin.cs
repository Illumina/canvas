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

        public GenomicBin(string chromosome, BedInterval interval)
        {
            Chromosome = chromosome;
            Interval = interval;
        }

        public GenomicBin(string chromosome, BedInterval interval, int gc)
        {
            Chromosome = chromosome;
            Interval = interval;
            GC = gc;
        }

        public string Chromosome { get; set; }
        public int GC { get; set; }
        public BedInterval Interval { get; set; }
    }

    /// <summary>
    /// Container that holds information about a genomic interval and how many reads were observed therein.
    /// </summary>
    public class SampleGenomicBin
    {
        private readonly GenomicBin _genomicBin;

        public float Count { get; set; }
        public double CountDeviation { get; set; }

        public SampleGenomicBin(string chr, int start, int stop, int gc, float count, double MadOfDIffs)
        {
            _genomicBin = new GenomicBin();
            GenomicBin.Chromosome = chr;
            GenomicBin.Interval = new BedInterval(start, stop);
            GenomicBin.GC = gc;
            this.Count = count;
            this.CountDeviation = MadOfDIffs;
        }

        public SampleGenomicBin(string chr, int start, int stop, int gc)
        {
            _genomicBin = new GenomicBin();
            GenomicBin.Chromosome = chr;
            GenomicBin.Interval = new BedInterval(start, stop);
            GenomicBin.GC = gc;
            this.CountDeviation = -1;
        }

        public SampleGenomicBin(string chr, int start, int stop, float count)
        {
            _genomicBin = new GenomicBin();
            GenomicBin.Chromosome = chr;
            GenomicBin.Interval = new BedInterval(start, stop);
            this.Count = count;
        }

        public SampleGenomicBin(string chr, int start, int stop, int gc, float count)
        {
            _genomicBin = new GenomicBin(chr, new BedInterval(start, stop), gc);
            this.Count = count;
            this.CountDeviation = -1;
        }

        public int Size
        {
            get { return GenomicBin.Interval.End - GenomicBin.Interval.Start; }
        }

        /// <summary>
        /// zero based inclusive start
        /// </summary>
        public int Start
        {
            get { return GenomicBin.Interval.Start; }
        }

        /// <summary>
        /// one-based inclusive end (or eqivalently the zero based exclusive end)
        /// </summary>
        public int Stop
        {
            get { return GenomicBin.Interval.End; }
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
