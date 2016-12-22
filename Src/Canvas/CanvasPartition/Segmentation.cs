using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using CanvasCommon;
using Isas.SequencingFiles;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Random;
using Utilities = CanvasCommon.Utilities;

namespace CanvasPartition
{
    class Segmentation
    {
        #region Members

        private string InputBinPath;
        private const int idxChr = 0, idxStart = 1, idxEnd = 2;
        private int idxScore = 3;
        public Dictionary<string, uint[]> StartByChr = new Dictionary<string, uint[]>();
        public Dictionary<string, uint[]> EndByChr = new Dictionary<string, uint[]>();
        public Dictionary<string, double[]> ScoreByChr = new Dictionary<string, double[]>();
        public string ForbiddenIntervalBedPath = null;
        public int MaxInterBinDistInSegment;
        #endregion

        public class Segment
        {
            public uint start; // Genomic start location
            public uint end; // Genomic end location (inclusive)
        }

        public class GenomeSegmentationResults
        {
            public IDictionary<string, Segmentation.Segment[]> SegmentByChr;

            public GenomeSegmentationResults(IDictionary<string, Segment[]> segmentByChr)
            {
                this.SegmentByChr = segmentByChr;
            }
        }

        public enum SegmentationMethod
        {
            Wavelets,
            CBS,
            HMM
        }

        public Segmentation(string inputBinPath, string forbiddenBedPath, int maxInterBinDistInSegment, 
            string dataType = "logratio")
                    {
            this.InputBinPath = inputBinPath;
            this.ForbiddenIntervalBedPath = forbiddenBedPath;
            this.MaxInterBinDistInSegment = maxInterBinDistInSegment;
            this.ReadBEDInput();
                    }

        private Segmentation.GenomeSegmentationResults GetDummySegmentationResults()
                    {
            Segmentation.GenomeSegmentationResults results = new Segmentation.GenomeSegmentationResults(new Dictionary<string, Segmentation.Segment[]>());
            return results;
            }


        public static Segmentation.Segment[] DeriveSegments(List<int> breakpoints, int sizeScoreByChr, uint[] startByChr,  uint[] endByChr)
            {
                    List<int> startBreakpointsPos = new List<int>();
                    List<int> endBreakpointPos = new List<int>();
                    List<int> lengthSeg = new List<int>();
                    if (breakpoints.Count() >= 2 && sizeScoreByChr > 10)
                    {
                        startBreakpointsPos.Add(breakpoints[0]);
                        endBreakpointPos.Add(breakpoints[1] - 1);
                        lengthSeg.Add(breakpoints[1] - 1);
                        for (int i = 1; i < breakpoints.Count - 1; i++)
                        {
                    startBreakpointsPos.Add(breakpoints[i]-1);
                    endBreakpointPos.Add(breakpoints[i + 1]);
                            lengthSeg.Add(breakpoints[i + 1] - 1 - breakpoints[i]);
                        }
                        startBreakpointsPos.Add(breakpoints[breakpoints.Count - 1]);
                        endBreakpointPos.Add(sizeScoreByChr - 1);
                        lengthSeg.Add(sizeScoreByChr - breakpoints[breakpoints.Count - 1] - 1);
                    }
                    else
                    {
                        startBreakpointsPos.Add(0);
                        endBreakpointPos.Add(sizeScoreByChr - 1);
                        lengthSeg.Add(sizeScoreByChr - 1);
                    }


            Segmentation.Segment[] segments = new Segmentation.Segment[startBreakpointsPos.Count];
                    for (int i = 0; i < startBreakpointsPos.Count; i++)
                    {
                        int start = startBreakpointsPos[i];
                        int end = endBreakpointPos[i];
                segments[i] = new Segmentation.Segment
                    {
                    start = startByChr[start],
                    end = endByChr[end]
                };
                // Genomic start
                // Genomic end
        }
            return segments;
                    }


        /// <summary>
        /// Assume that the rows are sorted by the start position and ascending order
        /// </summary>
        private void ReadBEDInput()
        {
            GenomicBinFilter binFilter = new GenomicBinFilter(ForbiddenIntervalBedPath);

            try
            {
                Dictionary<string, List<uint>> startByChr = new Dictionary<string, List<uint>>(),
                    endByChr = new Dictionary<string, List<uint>>();
                Dictionary<string, List<double>> scoreByChr = new Dictionary<string, List<double>>();
                using (GzipReader reader = new GzipReader(this.InputBinPath))
                {
                    string line;
                    string[] tokens;
                    while ((line = reader.ReadLine()) != null)
                    {
                        tokens = line.Split('\t');
                        string chrom = tokens[Segmentation.idxChr].Trim();
                        uint start = Convert.ToUInt32(tokens[Segmentation.idxStart].Trim());
                        uint end = Convert.ToUInt32(tokens[Segmentation.idxEnd].Trim());
                        if (binFilter.SkipBin(chrom, start, end))
                            continue;
                        if (!startByChr.ContainsKey(chrom))
                        {
                            startByChr.Add(chrom, new List<uint>());
                            endByChr.Add(chrom, new List<uint>());
                            scoreByChr.Add(chrom, new List<double>());
                        }
                        startByChr[chrom].Add(start);
                        endByChr[chrom].Add(end);
                        scoreByChr[chrom].Add(Convert.ToDouble(tokens[this.idxScore].Trim()));
                    }
                    foreach (string chr in startByChr.Keys)
                    {
                        this.StartByChr[chr] = startByChr[chr].ToArray();
                        this.EndByChr[chr] = endByChr[chr].ToArray();
                        this.ScoreByChr[chr] = scoreByChr[chr].ToArray();
                    }

                }
            }
            catch (Exception e)
            {
                Console.Error.WriteLine("File {0} could not be read:", this.InputBinPath);
                Console.Error.WriteLine(e.Message);
                Environment.Exit(1);
            }
        }


        /// <summary>
        /// Remap index from genomic coordinates into CanvasBin coordinates
        /// </summary>
        private static int? RemapIndex(uint[] startPos, uint[] endPos, int value, int length, ref int index)
        {
            const int distanceThreshold = 10000;
            var bestMinDistanceStart = Int32.MaxValue;
            var remappedIndex = 0;
            while (index < length)
            {
                int tmpMinDistanceStart = Math.Abs(Convert.ToInt32(startPos[index] + (endPos[index] - startPos[index])) - value);
                index++;
                if (tmpMinDistanceStart < bestMinDistanceStart)
                {
                    bestMinDistanceStart = tmpMinDistanceStart;
                    remappedIndex = index - 1;
                }
                else if (bestMinDistanceStart < distanceThreshold)
                {
                    return remappedIndex;
                }
                else
                {
                    return null;
                }
            }
            return null;
        }

        /// <summary>
        /// Remap GenomicBin from genome coordiantes into CanvasBin coordiantes
        /// </summary>
        public static List<SampleGenomicBin> RemapCommonRegions(List<SampleGenomicBin> commonRegions, uint[] startByChr, uint[] endByChr)
        {
            var length = startByChr.Length;
            var index = 0;
            List<SampleGenomicBin> commonRegionsRemapped = new List<SampleGenomicBin>();

            foreach (SampleGenomicBin commonRegion in commonRegions)
            {
                if (index > length)
                    break;
                var startSegment = RemapIndex(startByChr, endByChr, commonRegion.Start, length, ref index);
                var endSegment   = RemapIndex(startByChr, endByChr, commonRegion.Stop,  length, ref index);

                if (startSegment.HasValue && endSegment.HasValue)
                {
                    SampleGenomicBin interval = new SampleGenomicBin();
                    interval.Start = startSegment.Value;
                    interval.Stop = endSegment.Value;
                    commonRegionsRemapped.Add(interval);
                }
            }
            return commonRegionsRemapped;
        }

        /// <summary>
        /// Merge segmentation breakpoints with common CNV intervals
        /// </summary>
        public static List<int> OverlapCommonRegions(List<int> breakpoints, List<SampleGenomicBin> commonCNVintervals)
        {
            List<int> newBreakpoints = new List<int>();
            int index = 0;
            int length = commonCNVintervals.Count;
            foreach (int breakpoint in breakpoints)
            {
                while (index < length)
                {

                    if (breakpoint <= commonCNVintervals[index].Start)
                    {
                        newBreakpoints.Add(breakpoint);
                        break;
                    }
                    else if (breakpoint > commonCNVintervals[index].Start && breakpoint < commonCNVintervals[index].Stop)
                    {
                        newBreakpoints.Add(commonCNVintervals[index].Start);
                        newBreakpoints.Add(commonCNVintervals[index].Stop);
                        index++;
                        break;
                    }
                    else if (breakpoint >= commonCNVintervals[index].Stop)
                    {
                        newBreakpoints.Add(commonCNVintervals[index].Start);
                        newBreakpoints.Add(commonCNVintervals[index].Stop);
                        index++;
                    }
                }
                if (index > length)
                    newBreakpoints.Add(breakpoint);
            }
            return newBreakpoints;
        }


        public void WriteCanvasPartitionResults(string outPath, GenomeSegmentationResults regmentationResults)
        {
            Dictionary<string, bool> starts = new Dictionary<string, bool>();
            Dictionary<string, bool> stops = new Dictionary<string, bool>();

            foreach (string chr in regmentationResults.SegmentByChr.Keys)
            {
                for (int segmentIndex = 0; segmentIndex < regmentationResults.SegmentByChr[chr].Length; segmentIndex++)
                {
                    Segmentation.Segment segment = regmentationResults.SegmentByChr[chr][segmentIndex];
                    starts[chr + ":" + segment.start] = true;
                    stops[chr + ":" + segment.end] = true;
                }
            }

            Dictionary<string, List<SampleGenomicBin>> excludedIntervals = new Dictionary<string, List<SampleGenomicBin>>();
            if (!string.IsNullOrEmpty(ForbiddenIntervalBedPath))
            {
                excludedIntervals = CanvasCommon.Utilities.LoadBedFile(ForbiddenIntervalBedPath);
            }

            using (GzipWriter writer = new GzipWriter(outPath))
            {
                int segmentNum = -1;

                foreach (string chr in StartByChr.Keys)
                {
                    List<SampleGenomicBin> excludeIntervals = null;
                    if (excludedIntervals.ContainsKey(chr)) excludeIntervals = excludedIntervals[chr];
                    int excludeIndex = 0; // Points to the first interval which *doesn't* end before our current position
                    uint previousBinEnd = 0;
                    for (int pos = 0; pos < StartByChr[chr].Length; pos++)
                    {
                        uint start = StartByChr[chr][pos];
                        uint end = EndByChr[chr][pos];
                        string key = chr + ":" + start;
                        bool newSegment = IsNewSegment(starts, key, excludeIntervals, previousBinEnd, end, start, ref excludeIndex);
                        if (newSegment) segmentNum++;
                        writer.WriteLine(string.Format($"{chr}\t{start}\t{end}\t{ScoreByChr[chr][pos]}\t{segmentNum}"));
                        previousBinEnd = end;
                    }
                }
            }
        }

        private bool IsNewSegment(Dictionary<string, bool> starts, string key, List<SampleGenomicBin> excludeIntervals, uint previousBinEnd,
            uint end, uint start, ref int excludeIndex)
                        {
            bool newSegment = starts.ContainsKey(key) ? true : false;

                        if (excludeIntervals != null)
                        {
                while (excludeIndex < excludeIntervals.Count && excludeIntervals[excludeIndex].Stop < previousBinEnd)
                    excludeIndex++;
                            if (excludeIndex < excludeIntervals.Count)
                            {
                                // Note: forbiddenZoneMid should never fall inside a bin, becuase these intervals were already excluded 
                                // from consideration during the call to CanvasBin.
                    int forbiddenZoneMid = (excludeIntervals[excludeIndex].Start + excludeIntervals[excludeIndex].Stop)/2;
                                if (previousBinEnd < forbiddenZoneMid && end >= forbiddenZoneMid) newSegment = true;
                            }
                        }
                        if (previousBinEnd > 0 && MaxInterBinDistInSegment >= 0 && previousBinEnd + MaxInterBinDistInSegment < start
                            && !newSegment)
                        {
                            newSegment = true;
                        }
            return newSegment;
        }
    }
}
