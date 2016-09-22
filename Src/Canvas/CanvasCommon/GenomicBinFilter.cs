using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CanvasCommon
{
    public class GenomicBinFilter
    {
        private readonly Dictionary<string, List<GenomicBin>> _excludedIntervals;
        private string _prevChrom = null;
        private uint _prevStart = 0;
        private List<GenomicBin> _intervals = new List<GenomicBin>();
        private int _intervalIndex = -1;

        public GenomicBinFilter(string forbiddenIntervalBedPath)
        {
            _excludedIntervals = new Dictionary<string, List<GenomicBin>>();
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
                    new List<GenomicBin>();
                _intervalIndex = 0;
            }
            else if (start < _prevStart)
            {
                _intervalIndex = 0;
            }
            _prevStart = start;
            // while |- interval -|  |- bin -|
            while (_intervalIndex < _intervals.Count && _intervals[_intervalIndex].Stop < start + 1)
            {
                _intervalIndex++;
            }
            // now we either run out of intervals or  |- interval -| or           |- interval -|
            //                                            |- bin -|      |- bin -|
            // skip bins overlapping exclude intervals
            if (_intervalIndex < _intervals.Count && _intervals[_intervalIndex].Start < stop)
            {
                return true;
            }

            return false;
        }

        public bool SkipBin(GenomicBin bin)
        {
            return SkipBin(bin.Chromosome, (uint)bin.Start, (uint)bin.Stop);
        }
    }
}
