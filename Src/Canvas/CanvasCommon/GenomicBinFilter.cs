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
        private string _currChrom = null;
        private List<GenomicBin> _intervals = null;
        private int _intervalIndex = -1;

        public GenomicBinFilter(string forbiddenIntervalBedPath)
        {
            _excludedIntervals = new Dictionary<string, List<GenomicBin>>();
            if (!string.IsNullOrEmpty(forbiddenIntervalBedPath))
            {
                _excludedIntervals = Utilities.LoadBedFile(forbiddenIntervalBedPath);
            }
            Reset();
        }

        public void Reset()
        {
            _currChrom = null;
            _intervals = null;
            _intervalIndex = -1;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="chrom">chromosome</param>
        /// <param name="start">0-based start, inclusive</param>
        /// <param name="stop">0-based end, exclusive</param>
        /// <returns></returns>
        public bool SkipBin(string chrom, uint start, uint stop)
        {
            if (chrom != _currChrom)
            {
                _currChrom = chrom;
                _intervals = _excludedIntervals.ContainsKey(chrom) ? _excludedIntervals[chrom] : null;
                _intervalIndex = 0;
            }
            while (_intervals != null && _intervalIndex < _intervals.Count && _intervals[_intervalIndex].Stop < start + 1)
            {
                _intervalIndex++;
            }
            // skip bins overlapping exclude intervals
            if (_intervals != null && _intervalIndex < _intervals.Count && _intervals[_intervalIndex].Start < stop)
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
