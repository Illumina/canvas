using System;
using System.Collections.Generic;
using System.Linq;
using Isas.SequencingFiles;
using CanvasCommon;

namespace CanvasPartition
{
    internal class VafContainingBins
    {
        public int Index { get; }
        public double Vaf { get; }

        public VafContainingBins(int index, double vaf)
        {
            Index = index;
            Vaf = vaf;
        }
    }
    class SegmentationInput
    {
        public string ReferenceFolder { get; }

        #region Members

        private string InputBinPath;
        private string InputVafPath;
        private const int idxChr = 0, idxStart = 1, idxEnd = 2;
        private int idxScore = 3;
        public Dictionary<string, uint[]> StartByChr = new Dictionary<string, uint[]>();
        public Dictionary<string, uint[]> EndByChr = new Dictionary<string, uint[]>();
        public Dictionary<string, double[]> CoverageByChr = new Dictionary<string, double[]>();
        public Dictionary<string, List<VafContainingBins>> VafByChr = new Dictionary<string, List<VafContainingBins>>();

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
            public IDictionary<string, SegmentationInput.Segment[]> SegmentByChr;

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

        public SegmentationInput(string inputBinPath, string inputVafPath, string forbiddenBedPath, int maxInterBinDistInSegment,
            string referenceFolder, string dataType = "logratio")
        {
            InputBinPath = inputBinPath;
            InputVafPath = inputVafPath;
            ForbiddenIntervalBedPath = forbiddenBedPath;
            MaxInterBinDistInSegment = maxInterBinDistInSegment;
            ReadInputFiles(referenceFolder);
        }

        private void ReadInputFiles(string referenceFolder)
        {
            ReadBEDInput();
            LoadVAFInput(referenceFolder);
        }

        private SegmentationInput.GenomeSegmentationResults GetDummySegmentationResults()
        {
            var results = new SegmentationInput.GenomeSegmentationResults(new Dictionary<string, SegmentationInput.Segment[]>());
            return results;
        }


        public static SegmentationInput.Segment[] DeriveSegments(List<int> breakpoints, int segmentsLength, uint[] startByChr, uint[] endByChr)
        {
            if (segmentsLength < 1)
                return new Segment[0];
            var startBreakpointsPos = new List<int>();
            var endBreakpointPos = new List<int>();
            var lengthSeg = new List<int>();
            if (breakpoints.Count >= 2 && segmentsLength > 10)
            {
                startBreakpointsPos.Add(breakpoints[0]);
                endBreakpointPos.Add(breakpoints[1] - 1);
                lengthSeg.Add(breakpoints[1] - 1);
                for (int i = 1; i < breakpoints.Count - 1; i++)
                {
                    startBreakpointsPos.Add(breakpoints[i] - 1);
                    endBreakpointPos.Add(breakpoints[i + 1]);
                    lengthSeg.Add(breakpoints[i + 1] - 1 - breakpoints[i]);
                }
                startBreakpointsPos.Add(breakpoints[breakpoints.Count - 1]);
                endBreakpointPos.Add(segmentsLength - 1);
                lengthSeg.Add(segmentsLength - breakpoints[breakpoints.Count - 1] - 1);
            }
            else
            {
                startBreakpointsPos.Add(0);
                endBreakpointPos.Add(segmentsLength - 1);
                lengthSeg.Add(segmentsLength - 1);
            }


            var segments = new SegmentationInput.Segment[startBreakpointsPos.Count];
            for (int i = 0; i < startBreakpointsPos.Count; i++)
            {
                int start = startBreakpointsPos[i];
                int end = endBreakpointPos[i];
                segments[i] = new SegmentationInput.Segment
                {
                    start = startByChr[start],
                    end = endByChr[end]
                };
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
                var startByChr = new Dictionary<string, List<uint>>();
                var endByChr = new Dictionary<string, List<uint>>();
                var scoreByChr = new Dictionary<string, List<double>>();
                using (var reader = new GzipReader(this.InputBinPath))
                {
                    string line;
                    while ((line = reader.ReadLine()) != null)
                    {
                        var tokens = line.Split('\t');
                        string chrom = tokens[idxChr].Trim();
                        uint start = Convert.ToUInt32(tokens[idxStart].Trim());
                        uint end = Convert.ToUInt32(tokens[idxEnd].Trim());
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
                        scoreByChr[chrom].Add(Convert.ToDouble(tokens[idxScore].Trim()));
                    }
                    foreach (string chr in startByChr.Keys)
                    {
                        this.StartByChr[chr] = startByChr[chr].ToArray();
                        this.EndByChr[chr] = endByChr[chr].ToArray();
                        this.CoverageByChr[chr] = scoreByChr[chr].ToArray();
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
        /// Parse the outputs of CanvasSNV, and note these variant frequencies in the appropriate segment.
        /// </summary>
        public void LoadVAFInput(string referenceFolder)
        {
            try
            {
                var vafByChr = new Dictionary<string, List<List<double>>>();
                var intervalsByChromosome = new Dictionary<string, List<Interval>>();

                foreach (string chr in StartByChr.Keys)
                {
                    vafByChr[chr] = new List<List<double>>(StartByChr[chr].Length);
                    intervalsByChromosome[chr] = new List<Interval>();
                    for (int index = 0; index < StartByChr[chr].Length; index++)
                    {
                        vafByChr[chr].Add(new List<double>());
                        intervalsByChromosome[chr].Add(new Interval(StartByChr[chr][index], EndByChr[chr][index]));
                    }
                }

                var alleleCountsByChromosome = CanvasIO.ReadFrequencies(this.InputVafPath, intervalsByChromosome, referenceFolder, out float meanCoverage);

                foreach (var chr in alleleCountsByChromosome.Keys)
                {
                    for (int index = 0; index < alleleCountsByChromosome[chr].Count; index++)
                    {
                        vafByChr[chr][index] = alleleCountsByChromosome[chr][index].Select(genotype =>
                        Math.Max(genotype.CountsA, genotype.CountsB) / Convert.ToDouble(genotype.CountsA + genotype.CountsB)).ToList();
                    }
                }

                foreach (string chr in vafByChr.Keys)
                {
                    VafByChr[chr] = new List<VafContainingBins>();
                    var index = 0;
                    foreach (var bin in vafByChr[chr])
                    {
                        if (bin.Count < 1) continue;
                        VafByChr[chr].Add(new VafContainingBins(index, bin.Average()));
                        index++;
                    }
                }

            }
            catch (Exception e)
            {
                Console.Error.WriteLine("File {0} could not be read:", this.InputVafPath);
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
            int bestMinDistanceStart = Int32.MaxValue;
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
            int length = startByChr.Length;
            var index = 0;
            var commonRegionsRemapped = new List<SampleGenomicBin>();

            foreach (var commonRegion in commonRegions)
            {
                if (index > length)
                    break;
                var startSegmentIndex = RemapIndex(startByChr, endByChr, commonRegion.Start, length, ref index);
                var endSegmentIndex = RemapIndex(startByChr, endByChr, commonRegion.Stop, length, ref index);

                if (!startSegmentIndex.HasValue || !endSegmentIndex.HasValue) continue;

                var interval = new SampleGenomicBin
                {
                    Start = startSegmentIndex.Value,
                    Stop = endSegmentIndex.Value
                };
                commonRegionsRemapped.Add(interval);
            }
            return commonRegionsRemapped;
        }

        /// <summary>
        /// Merge segmentation breakpoints with common CNV intervals
        /// </summary>
        public static List<int> OverlapCommonRegions(List<int> breakpoints, List<SampleGenomicBin> commonCNVintervals)
        {
            var newBreakpoints = new List<int>();
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
                    if (breakpoint > commonCNVintervals[index].Start && breakpoint < commonCNVintervals[index].Stop)
                    {
                        newBreakpoints.Add(commonCNVintervals[index].Start);
                        newBreakpoints.Add(commonCNVintervals[index].Stop);
                        index++;
                        break;
                    }
                    if (breakpoint >= commonCNVintervals[index].Stop)
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


        public void WriteCanvasPartitionResults(string outPath, GenomeSegmentationResults segmentationResults)
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

            var excludedIntervals = new Dictionary<string, List<SampleGenomicBin>>();
            if (!string.IsNullOrEmpty(ForbiddenIntervalBedPath))
            {
                excludedIntervals = CanvasCommon.Utilities.LoadBedFile(ForbiddenIntervalBedPath);
            }

            using (var writer = new GzipWriter(outPath))
            {
                int segmentNum = -1;

                foreach (string chr in StartByChr.Keys)
                {
                    List<SampleGenomicBin> excludeIntervals = null;
                    if (excludedIntervals.ContainsKey(chr)) excludeIntervals = excludedIntervals[chr];
                    var excludeIndex = 0; // Points to the first interval which *doesn't* end before our current position
                    uint previousBinEnd = 0;
                    for (int pos = 0; pos < StartByChr[chr].Length; pos++)
                    {
                        uint start = StartByChr[chr][pos];
                        uint end = EndByChr[chr][pos];
                        string key = chr + ":" + start;
                        bool newSegment = IsNewSegment(starts, key, excludeIntervals, previousBinEnd, end, start, ref excludeIndex);
                        if (newSegment) segmentNum++;
                        writer.WriteLine(string.Format($"{chr}\t{start}\t{end}\t{CoverageByChr[chr][pos]}\t{segmentNum}"));
                        previousBinEnd = end;
                    }
                }
            }
        }

        private bool IsNewSegment(Dictionary<string, bool> starts, string key, List<SampleGenomicBin> excludeIntervals, uint previousBinEnd,
            uint end, uint start, ref int excludeIndex)
        {
            bool newSegment = starts.ContainsKey(key);

            if (excludeIntervals != null)
            {
                while (excludeIndex < excludeIntervals.Count && excludeIntervals[excludeIndex].Stop < previousBinEnd)
                    excludeIndex++;
                if (excludeIndex < excludeIntervals.Count)
                {
                    // Note: forbiddenZoneMid should never fall inside a bin, becuase these intervals were already excluded 
                    // from consideration during the call to CanvasBin.
                    int forbiddenZoneMid = (excludeIntervals[excludeIndex].Start + excludeIntervals[excludeIndex].Stop) / 2;
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

        private double GetEvennessScore()
        {
            throw new NotImplementedException();
        }
    }
}
