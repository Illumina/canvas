using System.Collections.Generic;

namespace CanvasCommon
{
    public class GenomicBinFilter
    {
        private readonly Dictionary<string, List<SampleGenomicBin>> _excludedIntervals;
        private string _prevChrom = null;
        private uint _prevStart = 0;
        private List<SampleGenomicBin> _intervals = new List<SampleGenomicBin>();
        private int _intervalIndex = -1;

        public GenomicBinFilter(string forbiddenIntervalBedPath)
        {
            _excludedIntervals = new Dictionary<string, List<SampleGenomicBin>>();
            if (!string.IsNullOrEmpty(forbiddenIntervalBedPath))
            {
                _excludedIntervals = Utilities.LoadBedFile(forbiddenIntervalBedPath);
            }
        }

        /// <summary>
        /// Assumes that the method is called on bins sorted by start in ascending order for each chromosome.
        /// </summary>
        /// <param name="chrom">chromosome</param>
        /// <param name="start">0-based start, inclusive</param>
        /// <param name="stop">0-based end, exclusive</param>
        /// <returns></returns>
        public bool SkipBin(string chrom, uint start, uint stop)
        {
            if (chrom != _prevChrom)
            {
                _prevChrom = chrom;
                _intervals = _excludedIntervals.ContainsKey(chrom) ? _excludedIntervals[chrom] :
                    new List<SampleGenomicBin>();
                _intervalIndex = 0;
            }
            else if (start < _prevStart)
            {
                _intervalIndex = 0;
            }
            _prevStart = start;

            for (; _intervalIndex < _intervals.Count; _intervalIndex++)
            {
                if (_intervals[_intervalIndex].Stop <= start) // |- interval -|  |- bin -|
                {
                    continue;
                }
                else if (_intervals[_intervalIndex].Start >= stop) // |- bin -|  |- interval -|
                {
                    return false;
                }
                return true; // overlaps
            }

            return false; // ran out of intervals
        }

        public bool SkipBin(SampleGenomicBin bin)
        {
            return SkipBin(bin.GenomicBin.Chromosome, (uint)bin.Start, (uint)bin.Stop);
        }
    }
}
