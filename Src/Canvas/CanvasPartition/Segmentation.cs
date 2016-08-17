using System;
using System.Collections.Generic;
using System.Collections.Concurrent;
using System.Runtime.InteropServices;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Diagnostics;
using SequencingFiles;
using MathNet.Numerics.Random;
using CanvasCommon;
using Isas.Shared;
using Utilities = CanvasCommon.Utilities;

namespace CanvasPartition
{
    class Segmentation
    {
        #region Members
        public static readonly double DefaultAlpha = 0.01; // default alpha for CBS
        public static readonly double DefaultMadFactor = 2.0; // default MAD factor for Wavelets
        private string InputBinPath;
        private string DataType;
        private const int idxChr = 0, idxStart = 1, idxEnd = 2;
        private int idxScore = 3;
        private Dictionary<string, uint[]> StartByChr = new Dictionary<string, uint[]>();
        private Dictionary<string, uint[]> EndByChr = new Dictionary<string, uint[]>();
        private Dictionary<string, double[]> ScoreByChr = new Dictionary<string, double[]>();
        private GenomeSegmentationResults SegmentationResults;
        public double Alpha = DefaultAlpha;
        public double MadFactor = DefaultMadFactor;
        public SegmentSplitUndo UndoMethod = SegmentSplitUndo.None;
        private string ForbiddenIntervalBedPath = null;
        private int MaxInterBinDistInSegment = 1000000;
        #endregion

        // dataType: "logratio" (aCGH, ROMA, etc.) or "binary" (LOH)
        public Segmentation(string inputBinPath, string forbiddenBedPath, string dataType = "logratio",
            int maxInterBinDistInSegment = 1000000)
        {
            this.InputBinPath = inputBinPath;
            this.DataType = dataType;
            this.SegmentationResults = null;
            this.ForbiddenIntervalBedPath = forbiddenBedPath;
            this.MaxInterBinDistInSegment = maxInterBinDistInSegment;
            // Read the input file:
            this.ReadBEDInput();
        }

        public void SegmentGenome(string outPath, SegmentationMethod method, bool isGermline, string commonCNVsbedPath)
        {
            switch (method)
            {
                case SegmentationMethod.Wavelets:
                default:// use Wavelets if CBS is not selected       
                    Console.WriteLine("{0} Running Wavelet Partitioning", DateTime.Now);
                    this.Wavelets(isGermline, commonCNVsbedPath, madFactor: MadFactor, verbose: 2);
                    break;
                case SegmentationMethod.CBS:
                    Console.WriteLine("{0} Running CBS Partitioning", DateTime.Now);
                    this.CBS(verbose: 2);
                    break;
            }
            Console.WriteLine("{0} Write CanvasPartition results:", DateTime.Now);
            this.WriteCanvasPartitionResults(outPath);
            Console.WriteLine("{0} CanvasPartition results written out", DateTime.Now);
        }

        private GenomeSegmentationResults GetDummySegmentationResults()
        {
            GenomeSegmentationResults results = new GenomeSegmentationResults(new Dictionary<string, Segment[]>());
            return results;
        }

        /// <summary>
        /// Wavelets: unbalanced HAAR wavelets segmentation 
        /// </summary>
        /// <param name="threshold">wavelets coefficient threshold</param>
        private void Wavelets(bool isGermline, string commonCNVs, double thresholdLower = 5, double thresholdUpper = 80, double madFactor = 2, int minSize = 10, int verbose = 1)
        {
            Dictionary<string, int[]> inaByChr = new Dictionary<string, int[]>();
            Dictionary<string, double[]> finiteScoresByChr = new Dictionary<string, double[]>();

            List<ThreadStart> tasks = new List<ThreadStart>();
            foreach (KeyValuePair<string, double[]> scoreByChrKVP in ScoreByChr)
            {
                tasks.Add(new ThreadStart(() =>
                {
                    string chr = scoreByChrKVP.Key;
                    int[] ina;
                    Helper.GetFiniteIndices(scoreByChrKVP.Value, out ina); // not NaN, -Inf, Inf

                    double[] scores;
                    if (ina.Length == scoreByChrKVP.Value.Length)
                    {
                        scores = scoreByChrKVP.Value;
                    }
                    else
                    {
                        Helper.ExtractValues<double>(scoreByChrKVP.Value, ina, out scores);
                    }

                    lock (finiteScoresByChr)
                    {
                        finiteScoresByChr[chr] = scores;
                        inaByChr[chr] = ina;
                    }

                }));
            }
            Isas.Shared.Utilities.DoWorkParallelThreads(tasks);
            // Quick sanity-check: If we don't have any segments, then return a dummy result.
            int n = 0;
            foreach (var list in finiteScoresByChr.Values)
            {
                n += list.Length;
            }
            if (n == 0)
            {
                this.SegmentationResults = this.GetDummySegmentationResults();
                return;
            }

            Dictionary<string, Segment[]> segmentByChr = new Dictionary<string, Segment[]>();


            // load common CNV segments
            Dictionary<string, List<GenomicBin>> commonCNVintervals = null;
            if (commonCNVs != null)
            {
                commonCNVintervals = Utilities.LoadBedFile(commonCNVs);
                Utilities.SortAndOverlapCheck(commonCNVintervals, commonCNVs);
            }

            tasks = new List<ThreadStart>();
            foreach (string chr in ScoreByChr.Keys)
            {
                tasks.Add(new ThreadStart(() =>
                {
                    int[] ina = inaByChr[chr];
                    List<int> breakpoints = new List<int>();
                    int sizeScoreByChr = this.ScoreByChr[chr].Length;
                    if (sizeScoreByChr > minSize)
                    {
                        WaveletSegmentation.HaarWavelets(this.ScoreByChr[chr].ToArray(), thresholdLower, thresholdUpper,
                            breakpoints, isGermline, madFactor: madFactor);
                    }

                    List<int> startBreakpointsPos = new List<int>();
                    List<int> endBreakpointPos = new List<int>();
                    List<int> lengthSeg = new List<int>();

                    if (commonCNVs != null)
                    {
                        if (commonCNVintervals.ContainsKey(chr))
                        {
                            List<GenomicBin> commonCNVintervalsByChr = commonCNVintervals[chr];
                            List<int> commonBreakpointsStart = new List<int>();
                            List<int> commonBreakpointsEnd = new List<int>();
                            MergeCommonRegions(commonCNVintervalsByChr, this.StartByChr[chr], this.EndByChr[chr],
                                ref commonBreakpointsStart, ref commonBreakpointsEnd);
                            List<int> oldbreakpoints = breakpoints;
                            breakpoints = OverlapCommonRegions(oldbreakpoints, commonBreakpointsStart, commonBreakpointsEnd);
                        }
                    }

                    if (breakpoints.Count() >= 2 && sizeScoreByChr > 10)
                    {
                        startBreakpointsPos.Add(breakpoints[0]);
                        endBreakpointPos.Add(breakpoints[1] - 1);
                        lengthSeg.Add(breakpoints[1] - 1);

                        for (int i = 1; i < breakpoints.Count - 1; i++)
                        {
                            startBreakpointsPos.Add(breakpoints[i]);
                            endBreakpointPos.Add(breakpoints[i + 1] - 1);
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
                    // estimate segment means 

                    double[] segmentMeans = new double[lengthSeg.Count()];
                    int ss = 0, ee = 0;
                    for (int i = 0; i < lengthSeg.Count(); i++)
                    {
                        ee += lengthSeg[i];
                        // Works even if weights == null
                        segmentMeans[i] = Helper.WeightedAverage(this.ScoreByChr[chr], null, iStart: ss, iEnd: ee);
                        ss = ee;
                    }

                    Segment[] segments = new Segment[startBreakpointsPos.Count];
                    for (int i = 0; i < startBreakpointsPos.Count; i++)
                    {
                        int start = startBreakpointsPos[i];
                        int end = endBreakpointPos[i];
                        segments[i] = new Segment();
                        segments[i].start = this.StartByChr[chr][start]; // Genomic start
                        segments[i].end = this.EndByChr[chr][end]; // Genomic end
                        segments[i].nMarkers = lengthSeg[i];
                        segments[i].mean = segmentMeans[i];
                    }

                    lock (segmentByChr)
                    {
                        segmentByChr[chr] = segments;
                    }
                }));

            }
            Console.WriteLine("{0} Launching wavelet tasks", DateTime.Now);
            Isas.Shared.Utilities.DoWorkParallelThreads(tasks);
            Console.WriteLine("{0} Completed wavelet tasks", DateTime.Now);
            this.SegmentationResults = new GenomeSegmentationResults(segmentByChr);
            Console.WriteLine("{0} Segmentation results complete", DateTime.Now);
        }


        /// <summary>
        /// CBS: circular binary segmentation porting the R function segment in DNAcopy
        /// </summary>
        /// <param name="alpha">Now in this.Alpha</param>
        /// <param name="nPerm"></param>
        /// <param name="pMethod">"hybrid" or "perm"</param>
        /// <param name="minWidth"></param>
        /// <param name="kMax"></param>
        /// <param name="nMin"></param>
        /// <param name="eta"></param>
        /// <param name="sbdry"></param>
        /// <param name="trim"></param>
        /// <param name="undoSplit">"none" or "prune" or "sdundo"; now in this.UndoMethod</param>
        /// <param name="undoPrune"></param>
        /// <param name="undoSD"></param>
        /// <param name="verbose"></param>
        private void CBS(uint nPerm = 10000, string pMethod = "hybrid", int minWidth = 2, int kMax = 25,
            uint nMin = 200, double eta = 0.05, uint[] sbdry = null, double trim = 0.025,
            double undoPrune = 0.05, double undoSD = 3, int verbose = 1)
        {
            if (minWidth < 2 || minWidth > 5)
            {
                Console.Error.WriteLine("Minimum segment width should be between 2 and 5");
                Environment.Exit(1);
            }
            if (nMin < 4 * kMax)
            {
                Console.Error.WriteLine("nMin should be >= 4 * kMax");
                Environment.Exit(1);
            }
            if (sbdry == null)
            {
                GetBoundary.ComputeBoundary(nPerm, this.Alpha, eta, out sbdry);
            }

            Dictionary<string, int[]> inaByChr = new Dictionary<string, int[]>();
            Dictionary<string, double[]> finiteScoresByChr = new Dictionary<string, double[]>();

            List<ThreadStart> tasks = new List<ThreadStart>();
            foreach (KeyValuePair<string, double[]> scoreByChrKVP in ScoreByChr)
            {
                tasks.Add(new ThreadStart(() =>
                {
                    string chr = scoreByChrKVP.Key;
                    int[] ina;
                    Helper.GetFiniteIndices(scoreByChrKVP.Value, out ina); // not NaN, -Inf, Inf

                    double[] scores;
                    if (ina.Length == scoreByChrKVP.Value.Length)
                    {
                        scores = scoreByChrKVP.Value;
                    }
                    else
                    {
                        Helper.ExtractValues<double>(scoreByChrKVP.Value, ina, out scores);
                    }

                    lock (finiteScoresByChr)
                    {
                        finiteScoresByChr[chr] = scores;
                        inaByChr[chr] = ina;
                    }

                }));
            }
            Isas.Shared.Utilities.DoWorkParallelThreads(tasks);

            // Quick sanity-check: If we don't have any segments, then return a dummy result.
            int n = 0;
            foreach (var list in finiteScoresByChr.Values)
            {
                n += list.Length;
            }
            if (n == 0)
            {
                this.SegmentationResults = this.GetDummySegmentationResults();
                return;
            }

            double trimmedSD = Math.Sqrt(ChangePoint.TrimmedVariance(finiteScoresByChr, trim: trim));
            Dictionary<string, Segment[]> segmentByChr = new Dictionary<string, Segment[]>();

            // when parallelizing we need an RNG for each chromosome to get deterministic results
            Random seedGenerator = new MersenneTwister(0);
            Dictionary<string, Random> perChromosomeRandom = new Dictionary<string, Random>();
            foreach (string chr in this.ScoreByChr.Keys)
            {
                perChromosomeRandom[chr] = new MersenneTwister(seedGenerator.NextFullRangeInt32(), true);
            }

            tasks = new List<ThreadStart>();
            foreach (string chr in ScoreByChr.Keys)
            {
                tasks.Add(new ThreadStart(() =>
                {
                    int[] ina = inaByChr[chr];
                    int[] lengthSeg;
                    double[] segmentMeans;
                    ChangePoint.ChangePoints(this.ScoreByChr[chr], sbdry, out lengthSeg, out segmentMeans, perChromosomeRandom[chr],
                        dataType: this.DataType, alpha: this.Alpha, nPerm: nPerm,
                        pMethod: pMethod, minWidth: minWidth, kMax: kMax, nMin: nMin, trimmedSD: trimmedSD,
                        undoSplits: this.UndoMethod, undoPrune: undoPrune, undoSD: undoSD, verbose: verbose);

                    Segment[] segments = new Segment[lengthSeg.Length];
                    int cs1 = 0, cs2 = -1; // cumulative sum
                    for (int i = 0; i < lengthSeg.Length; i++)
                    {
                        cs2 += lengthSeg[i];
                        int start = ina[cs1];
                        int end = ina[cs2];
                        segments[i] = new Segment();
                        segments[i].start = this.StartByChr[chr][start]; // Genomic start
                        segments[i].end = this.EndByChr[chr][end]; // Genomic end
                        segments[i].nMarkers = lengthSeg[i];
                        segments[i].mean = segmentMeans[i];
                        cs1 += lengthSeg[i];
                    }

                    lock (segmentByChr)
                    {
                        segmentByChr[chr] = segments;
                    }
                }));
            }

            //Parallel.ForEach(tasks, t => { t.Invoke(); });
            Isas.Shared.Utilities.DoWorkParallelThreads(tasks);
            this.SegmentationResults = new GenomeSegmentationResults(segmentByChr);
        }

        private sealed class Segment
        {
            public uint start; // Genomic start location
            public uint end; // Genomic end location (inclusive)
            public int nMarkers; // Number of markers in a segment
            public double mean; // Mean signal acroos markers in a segment
        }

        private sealed class GenomeSegmentationResults
        {
            public IDictionary<string, Segment[]> SegmentByChr;

            public GenomeSegmentationResults(IDictionary<string, Segment[]> segmentByChr)
            {
                this.SegmentByChr = segmentByChr;
            }
        }

        /// <summary>
        /// Assume that the rows are sorted by the start position and ascending order
        /// </summary>
        private void ReadBEDInput()
        {
            try
            {
                Dictionary<string, List<uint>> startByChr = new Dictionary<string, List<uint>>(),
                    endByChr = new Dictionary<string, List<uint>>();
                Dictionary<string, List<double>> scoreByChr = new Dictionary<string, List<double>>();
                // Create an instance of StreamReader to read from a file. 
                // The using statement also closes the StreamReader. 
                using (GzipReader reader = new GzipReader(this.InputBinPath))
                {
                    string line;
                    string[] tokens;
                    while ((line = reader.ReadLine()) != null)
                    {
                        tokens = line.Split('\t');
                        string chr = tokens[Segmentation.idxChr].Trim();
                        if (!startByChr.ContainsKey(chr))
                        {
                            startByChr.Add(chr, new List<uint>());
                            endByChr.Add(chr, new List<uint>());
                            scoreByChr.Add(chr, new List<double>());
                        }
                        startByChr[chr].Add(Convert.ToUInt32(tokens[Segmentation.idxStart].Trim()));
                        endByChr[chr].Add(Convert.ToUInt32(tokens[Segmentation.idxEnd].Trim()));
                        scoreByChr[chr].Add(Convert.ToDouble(tokens[this.idxScore].Trim()));
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


        private void MergeCommonRegions(List<GenomicBin> commonRegions, uint[] startByChr, uint[] endByChr, ref List<int> bestMinDistanceStarts, ref List<int> bestMinDistanceStops)
        {
            int length = startByChr.Length;
            int indexByChr = 0;
            const int distanceThreshold = 10000;

            foreach (GenomicBin commonRegion in commonRegions)
            {

                if (indexByChr > length)
                    break;
                var bestMinDistanceStart = Int32.MaxValue;
                var startSegment = 0;
                while (indexByChr < length)
                {
                    int tmpMinDistanceStart = Math.Abs(Convert.ToInt32(startByChr[indexByChr] + (endByChr[indexByChr] - startByChr[indexByChr])) - commonRegion.Start);
                    indexByChr++;
                    if (tmpMinDistanceStart < bestMinDistanceStart)
                    {
                        bestMinDistanceStart = tmpMinDistanceStart;
                        startSegment = indexByChr - 1;
                    }
                    else
                        break;
                }
                var bestMinDistanceStop = Int32.MaxValue;
                var endSegment = 0;
                while (indexByChr < length)
                {
                    int tmpMinDistanceStop = Math.Abs(Convert.ToInt32(startByChr[indexByChr] + (endByChr[indexByChr] - startByChr[indexByChr])) - commonRegion.Stop);
                    indexByChr++;
                    if (tmpMinDistanceStop < bestMinDistanceStop)
                    {
                        bestMinDistanceStop = tmpMinDistanceStop;
                        endSegment = indexByChr - 1;
                    }               
                    else
                        break;
                }
                if (bestMinDistanceStart < distanceThreshold && bestMinDistanceStop < distanceThreshold)
                {
                    bestMinDistanceStarts.Add(startSegment);
                    bestMinDistanceStops.Add(endSegment);
                }
            }
        }

        public static List<int> OverlapCommonRegions(List<int> breakpoints, List<int> commonBreakpointsStart,
            List<int> commonBreakpointsEnd)
        {
            List<int> newBreakpoints = new List<int>();
            int commonBreakpointsIndex = 0;
            int length = commonBreakpointsStart.Count;
            foreach (int bkpt in breakpoints)
            {
                while (commonBreakpointsIndex < length)
                {

                    if(bkpt <= commonBreakpointsStart[commonBreakpointsIndex])
                    {
                        newBreakpoints.Add(bkpt);
                        break;
                    }
                    else if (bkpt > commonBreakpointsStart[commonBreakpointsIndex] && bkpt < commonBreakpointsEnd[commonBreakpointsIndex])
                    {
                        newBreakpoints.Add(commonBreakpointsStart[commonBreakpointsIndex]);
                        newBreakpoints.Add(commonBreakpointsEnd[commonBreakpointsIndex]);
                        commonBreakpointsIndex++;
                        break;
                    }
                    else if (bkpt >= commonBreakpointsEnd[commonBreakpointsIndex])
                    {
                        newBreakpoints.Add(commonBreakpointsStart[commonBreakpointsIndex]);
                        newBreakpoints.Add(commonBreakpointsEnd[commonBreakpointsIndex]);
                        commonBreakpointsIndex++;
                    }
                }
                if (commonBreakpointsIndex > length)
                    newBreakpoints.Add(bkpt);
            }
            return newBreakpoints;
        }

        private
            void WriteCanvasPartitionResults(string outPath)
        {
            Dictionary<string, bool> starts = new Dictionary<string, bool>();
            Dictionary<string, bool> stops = new Dictionary<string, bool>();

            foreach (string chr in SegmentationResults.SegmentByChr.Keys)
            {
                for (int segmentIndex = 0; segmentIndex < SegmentationResults.SegmentByChr[chr].Length; segmentIndex++)
                {
                    Segment segment = SegmentationResults.SegmentByChr[chr][segmentIndex];
                    starts[chr + ":" + segment.start] = true;
                    stops[chr + ":" + segment.end] = true;
                }
            }

            Dictionary<string, List<GenomicBin>> ExcludedIntervals = new Dictionary<string, List<GenomicBin>>();
            if (!string.IsNullOrEmpty(ForbiddenIntervalBedPath))
            {
                ExcludedIntervals = CanvasCommon.Utilities.LoadBedFile(ForbiddenIntervalBedPath);
            }

            using (GzipWriter writer = new GzipWriter(outPath))
            {
                int segmentNum = -1;

                foreach (string chr in StartByChr.Keys)
                {
                    List<GenomicBin> excludeIntervals = null;
                    if (ExcludedIntervals.ContainsKey(chr)) excludeIntervals = ExcludedIntervals[chr];
                    int excludeIndex = 0; // Points to the first interval which *doesn't* end before our current position
                    uint previousBinEnd = 0;
                    for (int pos = 0; pos < StartByChr[chr].Length; pos++)
                    {
                        uint start = StartByChr[chr][pos];
                        uint end = EndByChr[chr][pos];
                        bool newSegment = false;
                        string key = chr + ":" + start;
                        if (starts.ContainsKey(key))
                        {
                            newSegment = true;
                        }

                        if (excludeIntervals != null)
                        {
                            while (excludeIndex < excludeIntervals.Count && excludeIntervals[excludeIndex].Stop < previousBinEnd) excludeIndex++;
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

                        if (newSegment) segmentNum++;
                        writer.WriteLine(string.Format("{0}\t{1}\t{2}\t{3}\t{4}", chr, start, end, ScoreByChr[chr][pos], segmentNum));
                        previousBinEnd = end;
                    }
                }
            }
        }
    }
}
