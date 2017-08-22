using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using Isas.SequencingFiles;
using CanvasCommon;
using Illumina.Common.FileSystem;
using Isas.Framework.Logging;

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
        public string CoverageMetricsFile { get; }
        public string ReferenceFolder { get; }


        #region Members
        private string InputBinPath;
        private string InputVafPath;
        public Dictionary<string, List<VafContainingBins>> VafByChr = new Dictionary<string, List<VafContainingBins>>();
        public string ForbiddenIntervalBedPath = null;
        public int MaxInterBinDistInSegment;
        public CoverageInfo CoverageInfo = new CoverageInfo();
        ILogger _logger;
        #endregion

        public class Segment
        {
            public uint start; // Genomic start location
            public uint end; // Genomic end location (inclusive)
        }

        public class GenomeSegmentationResults
        {
            public IDictionary<string, SegmentationInput.Segment[]> SegmentByChr;

            public GenomeSegmentationResults(IDictionary<string, SegmentationInput.Segment[]> segmentByChr)
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
            string referenceFolder, string coverageMetricsFile, string dataType = "logratio")
        {
            CoverageMetricsFile = coverageMetricsFile;
            InputBinPath = inputBinPath;
            InputVafPath = inputVafPath;
            ForbiddenIntervalBedPath = forbiddenBedPath;
            MaxInterBinDistInSegment = maxInterBinDistInSegment;
            ReadInputFiles(referenceFolder);
        }

        private void ReadInputFiles(string referenceFolder)
        {
            CoverageInfo = CanvasSegment.ReadBEDInput(InputBinPath, ForbiddenIntervalBedPath);
            if (InputVafPath != null)
                LoadVAFInput(referenceFolder);
        }

        public static SegmentationInput.Segment[] DeriveSegments(List<int> breakpoints, int segmentsLength, uint[] startByChr, uint[] endByChr)
        {
            var startBreakpointsPos = new List<int>();
            var endBreakpointPos = new List<int>();
            var lengthSeg = new List<int>();
            if (breakpoints.Count >= 2 && segmentsLength > 10)
            {
                if (breakpoints[0] != 0)
                    breakpoints.Insert(0,0);
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
        /// Parse the outputs of CanvasSNV, and note these variant frequencies in the appropriate segment.
        /// </summary>
        public void LoadVAFInput(string referenceFolder)
        {
            try
            {
                var vafByChr = new Dictionary<string, List<List<double>>>();
                var intervalsByChromosome = new Dictionary<string, List<BedInterval>>();

                foreach (string chr in CoverageInfo.StartByChr.Keys)
                {
                    vafByChr[chr] = new List<List<double>>(CoverageInfo.StartByChr[chr].Length);
                    intervalsByChromosome[chr] = new List<BedInterval>();
                    for (int index = 0; index < CoverageInfo.StartByChr[chr].Length; index++)
                    {
                        vafByChr[chr].Add(new List<double>());
                        intervalsByChromosome[chr].Add(new BedInterval(Convert.ToInt32(CoverageInfo.StartByChr[chr][index]), 
                            Convert.ToInt32(CoverageInfo.EndByChr[chr][index])));
                    }
                }

                var alleleCountsByChromosome = CanvasIO.ReadFrequenciesWrapper(_logger, new FileLocation(this.InputVafPath), intervalsByChromosome);
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
                        if (bin.Count > 0) 
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
        /// Merge segmentation breakpoints with common CNV intervals
        /// </summary>
        public static List<int> OverlapCommonRegions(List<int> breakpoints, List<BedInterval> commonCNVintervals)
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
                    if (breakpoint > commonCNVintervals[index].Start && breakpoint < commonCNVintervals[index].End)
                    {
                        newBreakpoints.Add(commonCNVintervals[index].Start);
                        newBreakpoints.Add(commonCNVintervals[index].End);
                        index++;
                        break;
                    }
                    if (breakpoint >= commonCNVintervals[index].End)
                    {
                        newBreakpoints.Add(commonCNVintervals[index].Start);
                        newBreakpoints.Add(commonCNVintervals[index].End);
                        index++;
                    }
                }
                if (index > length)
                    newBreakpoints.Add(breakpoint);
            }
            return newBreakpoints;
        }


        public void WriteCanvasPartitionResults(string outPath, SegmentationInput.GenomeSegmentationResults segmentationResults)
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

                foreach (string chr in CoverageInfo.StartByChr.Keys)
                {
                    List<SampleGenomicBin> excludeIntervals = null;
                    if (excludedIntervals.ContainsKey(chr)) excludeIntervals = excludedIntervals[chr];
                    var excludeIndex = 0; // Points to the first interval which *doesn't* end before our current position
                    uint previousBinEnd = 0;
                    for (int pos = 0; pos < CoverageInfo.StartByChr[chr].Length; pos++)
                    {
                        uint start = CoverageInfo.StartByChr[chr][pos];
                        uint end = CoverageInfo.EndByChr[chr][pos];
                        string key = chr + ":" + start;
                        bool newSegment = IsNewSegment(starts, key, excludeIntervals, previousBinEnd, end, start, ref excludeIndex);
                        if (newSegment) segmentNum++;
                        writer.WriteLine(string.Format($"{chr}\t{start}\t{end}\t{CoverageInfo.CoverageByChr[chr][pos]}\t{segmentNum}"));
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

        /// <summary>
        /// Implements evenness score from https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkq072#55451628
        /// </summary>
        /// <returns></returns>
        public double GetEvennessScore(int windowSize)
        {
            const double IQRthreshold = 0.015;
            const int windowSizeIQR = 10000;
            var evennessScoresIQR = reportScoresByWindow(windowSizeIQR);
            var quartiles = CanvasCommon.Utilities.Quartiles(evennessScoresIQR.Select(Convert.ToSingle).ToList());
            var evennessScores = reportScoresByWindow(windowSize);
            double median = CanvasCommon.Utilities.Median(evennessScores.ToList());
            return quartiles.Item3 - quartiles.Item1 > IQRthreshold ? quartiles.Item3 * 100.0 : median * 100.0;
        }

        private ConcurrentBag<double> reportScoresByWindow(int windowSize)
        {
            var evennessScores = new ConcurrentBag<double>();
            var tasks = CoverageInfo.CoverageByChr.Select(coverage => new ThreadStart(() =>
            {
                for (var index = 0; index < coverage.Value.Length - windowSize; index += windowSize)
                {
                    var tmp = coverage.Value.Skip(index).Take(windowSize - 1).ToList();
                    double average = tmp.Average();
                    double tmpEvenness = 0;
                    for (var coverageBin = 0; coverageBin <= average; coverageBin++)
                    {
                        tmpEvenness += tmp.Select(bin => bin >= coverageBin).Count(c => c) / tmp.Sum();
                    }

                    if (!double.IsInfinity(tmpEvenness) && !double.IsNaN(tmpEvenness))
                        evennessScores.Add(tmpEvenness);
                }
            })).ToList();
            Parallel.ForEach(tasks, task => task.Invoke());
            return evennessScores;
        }
    }
}
