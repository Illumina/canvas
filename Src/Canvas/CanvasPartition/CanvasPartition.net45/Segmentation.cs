using System;
using System.Collections.Generic;
using System.Linq;
using Isas.SequencingFiles;
using CanvasCommon;

namespace CanvasPartition
{
    internal class CoverageToVafMapper
    {
        public int Index { get; }
        public double Vaf { get; }

        public CoverageToVafMapper(int index, double vaf)
        {
            Index = index;
            Vaf = vaf;
        }
    }
    class Segmentation
    {
        #region Members

        private string InputBinPath;
        private string InputVafPath;
        private const int idxChr = 0, idxStart = 1, idxEnd = 2;
        private int idxScore = 3;
        public Dictionary<string, uint[]> StartByChr = new Dictionary<string, uint[]>();
        public Dictionary<string, uint[]> EndByChr = new Dictionary<string, uint[]>();
        public Dictionary<string, double[]> CoverageByChr = new Dictionary<string, double[]>();
        public Dictionary<string, List<CoverageToVafMapper>> VafByChr = new Dictionary<string, List<CoverageToVafMapper>>();

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

        public Segmentation(string inputBinPath, string inputVafPath, string forbiddenBedPath, int maxInterBinDistInSegment,
            string dataType = "logratio")
        {
            InputBinPath = inputBinPath;
            InputVafPath = inputVafPath;
            ForbiddenIntervalBedPath = forbiddenBedPath;
            MaxInterBinDistInSegment = maxInterBinDistInSegment;
            ReadInputFiles();
        }

        private void ReadInputFiles()
        {
            ReadBEDInput();
            ReadVAFInput();
        }

        private Segmentation.GenomeSegmentationResults GetDummySegmentationResults()
        {
            var results = new Segmentation.GenomeSegmentationResults(new Dictionary<string, Segmentation.Segment[]>());
            return results;
        }


        public static Segmentation.Segment[] DeriveSegments(List<int> breakpoints, int sizeScoreByChr, uint[] startByChr, uint[] endByChr)
        {
            var startBreakpointsPos = new List<int>();
            var endBreakpointPos = new List<int>();
            var lengthSeg = new List<int>();
            if (breakpoints.Count >= 2 && sizeScoreByChr > 10)
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
                endBreakpointPos.Add(sizeScoreByChr - 1);
                lengthSeg.Add(sizeScoreByChr - breakpoints[breakpoints.Count - 1] - 1);
            }
            else
            {
                startBreakpointsPos.Add(0);
                endBreakpointPos.Add(sizeScoreByChr - 1);
                lengthSeg.Add(sizeScoreByChr - 1);
            }


            var segments = new Segmentation.Segment[startBreakpointsPos.Count];
            for (int i = 0; i < startBreakpointsPos.Count; i++)
            {
                int start = startBreakpointsPos[i];
                int end = endBreakpointPos[i];
                segments[i] = new Segmentation.Segment
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
        public void ReadVAFInput()
        {
            try
            {
                var vfByChr = new Dictionary<string, List<List<double>>>();
                foreach (string chr in StartByChr.Keys)
                {
                    vfByChr[chr] = new List<List<double>>(StartByChr[chr].Length);
                }
                Console.WriteLine("{0} Load variant frequencies from {1}", DateTime.Now, InputVafPath);

                using (var reader = new GzipReader(InputVafPath))
                {
                    while (true)
                    {
                        string fileLine = reader.ReadLine();
                        if (fileLine == null) break;
                        if (fileLine.Length == 0 || fileLine[0] == '#') continue; // Skip headers
                        var bits = fileLine.Split('\t');
                        if (bits.Length < 6)
                        {
                            Console.Error.WriteLine("* Bad line in {0}: '{1}'", InputVafPath, fileLine);
                            continue;
                        }
                        string chromosome = bits[0];
                        int position = int.Parse(bits[1]); // 1-based (from the input VCF to Canvas SNV)
                        int countRef = int.Parse(bits[4]);
                        int countAlt = int.Parse(bits[5]);
                        if (countRef + countAlt < 10) continue;
                        double VF = Math.Max(countRef, countAlt) / (double) (countRef + countAlt);
                        // Binary search for the segment this variant hits:
                        var start = 0;
                        int end = EndByChr[chromosome].Length - 1;
                        int mid = (start + end) / 2;
                        while (start <= end)
                        {
                            if (EndByChr[chromosome][mid] < position)
                            {
                                start = mid + 1;
                                mid = (start + end) / 2;
                                continue;
                            }
                            if (StartByChr[chromosome][mid] + 1 > position)
                            {
                                end = mid - 1;
                                mid = (start + end) / 2;
                                continue;
                            }
                            vfByChr[chromosome][mid].Add(VF);
                            break;
                        }
                    }
                    foreach (string chr in vfByChr.Keys)
                    {
                        VafByChr[chr] = vfByChr[chr].Where(bin => bin.Count > 1).
                            Select((bin, index) => new CoverageToVafMapper(index, bin.Average())).ToList();
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
                var startSegment = RemapIndex(startByChr, endByChr, commonRegion.Start, length, ref index);
                var endSegment = RemapIndex(startByChr, endByChr, commonRegion.Stop, length, ref index);

                if (!startSegment.HasValue || !endSegment.HasValue) continue;

                var interval = new SampleGenomicBin
                {
                    Start = startSegment.Value,
                    Stop = endSegment.Value
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
