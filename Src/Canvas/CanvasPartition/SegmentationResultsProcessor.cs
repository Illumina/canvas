using System.Collections.Generic;
using System.Diagnostics;
using CanvasCommon;
using CanvasPartition.Models;
using Isas.SequencingFiles;

namespace CanvasPartition
{
    internal class SegmentationResultsProcessor : ISegmentationResultsProcessor
    {
        private readonly int _maxInterBinDistInSegment;
        public SegmentationResultsProcessor(int maxInterBinDistInSegment)
        {
            _maxInterBinDistInSegment = maxInterBinDistInSegment;
        }

        public Dictionary<string, List<SegmentWithBins>> PostProcessSegments(
            GenomeSegmentationResults segmentationResults,
            PloidyInfo referencePloidy, Dictionary<string, List<SampleGenomicBin>> excludedIntervals, CoverageInfo coverageInfo)
        {
            var starts = new Dictionary<string, bool>();
            var stops = new Dictionary<string, bool>();

            foreach (string chr in segmentationResults.SegmentByChr.Keys)
            {
                for (int segmentIndex = 0; segmentIndex < segmentationResults.SegmentByChr[chr].Length; segmentIndex++)
                {
                    var segment = segmentationResults.SegmentByChr[chr][segmentIndex];
                    starts[chr + ":" + segment.start] = true;
                    stops[chr + ":" + segment.end] = true;
                }
            }

            int segmentNum = -1;


            var segmentsByChromosome = new Dictionary<string, List<SegmentWithBins>>();
            foreach (string chr in coverageInfo.StartByChr.Keys)
            {
                segmentsByChromosome.Add(chr, new List<SegmentWithBins>());
                SegmentWithBins currentSegment = null;

                List<SampleGenomicBin> excludeIntervals = null;

                if (excludedIntervals.ContainsKey(chr)) excludeIntervals = excludedIntervals[chr];
                var excludeIndex = 0; // Points to the first interval which *doesn't* end before our current position
                uint previousBinEnd = 0;

                for (int binIndex = 0; binIndex < coverageInfo.StartByChr[chr].Length; binIndex++)
                {
                    uint start = coverageInfo.StartByChr[chr][binIndex];
                    uint end = coverageInfo.EndByChr[chr][binIndex];

                    bool newSegment = IsNewSegment(starts, chr, excludeIntervals, previousBinEnd, end, start, ref excludeIndex, referencePloidy);

                    var bin = new Bin(start, end, coverageInfo.CoverageByChr[chr][binIndex]);
                    if (newSegment)
                    {
                        segmentNum++;
                        currentSegment = new SegmentWithBins(segmentNum, bin);
                        segmentsByChromosome[chr].Add(currentSegment);
                    }
                    else
                    {
                        if (currentSegment == null)
                        {
                            currentSegment = new SegmentWithBins(segmentNum, bin);
                            segmentsByChromosome[chr].Add(currentSegment);
                        }
                        else
                        {
                            currentSegment.AddBin(bin);
                        }
                    }


                    previousBinEnd = end;
                }

            }

            return segmentsByChromosome;

        }


        private bool IsNewSegment(Dictionary<string, bool> starts, string chr, List<SampleGenomicBin> excludeIntervals, uint previousBinEnd,
            uint end, uint start, ref int excludeIndex, PloidyInfo referencePloidy)
        {
            string key = chr + ":" + start;
            bool newSegment = starts.ContainsKey(key);

            if (excludeIntervals != null)
            {
                while (excludeIndex < excludeIntervals.Count && excludeIntervals[excludeIndex].Stop < previousBinEnd)
                    excludeIndex++;
                if (excludeIndex < excludeIntervals.Count)
                {
                    // Note: forbiddenZoneMid should never fall inside a bin, becuase these intervals were already excluded 
                    // from consideration during the call to CanvasBin.
                    var excludeStart = excludeIntervals[excludeIndex].Start;
                    var excludeStop = excludeIntervals[excludeIndex].Stop;
                    int forbiddenZoneMid = (excludeStart + excludeStop) / 2;
                    if (previousBinEnd < forbiddenZoneMid && end >= forbiddenZoneMid)
                    {
                        //Debug.Assert(previousBinEnd <= excludeStart);
                        //Debug.Assert(start >= excludeStop);
                        newSegment = true;
                    }
                }
            }
            if (previousBinEnd > 0 && _maxInterBinDistInSegment >= 0 && previousBinEnd + _maxInterBinDistInSegment < start
                && !newSegment)
            {
                newSegment = true;
            }
            // also start new segment if reference ploidy changes between end of last and end of this;
            // note that Interval takes 1-based positions, so using "previousBinEnd" effectively
            // includes the last base of the previous bin, allowing for a change at the bin boundary
            if (!newSegment && referencePloidy != null)
            {
                var refIval = new ReferenceInterval(chr, new Interval(previousBinEnd > 0 ? previousBinEnd : 1, end));
                if (!referencePloidy.IsUniformReferencePloidy(refIval))
                {
                    newSegment = true;
                }
            }
            return newSegment;
        }

    }
}