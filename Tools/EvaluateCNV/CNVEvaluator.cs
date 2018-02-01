using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Net.Sockets;
using CanvasCommon;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Illumina.Common.MathUtilities;
using Isas.Framework.Common_Statistics;

namespace EvaluateCNV
{
    public class BaseCounter
    {
        public int TotalVariants { get; set; }
        public long TotalVariantBases { get; set; }
        public double MeanAccuracy { get; set; }
        public double MedianAccuracy { get; set; }
        public int MinSize { get; }
        public int MaxSize { get; }
        // 3D array stores knownCn, CN call and REF ploidy
        public long[,,] BaseCount;
        public long[,,] RoiBaseCount;


        public BaseCounter(int maxCn, int minSize, int maxSize, bool hasRoi = false)
        {
            MinSize = minSize;
            MaxSize = maxSize;
            BaseCount = new long[maxCn + 1, maxCn + 1, 3];
            if (hasRoi)
                RoiBaseCount = new long[maxCn + 1, maxCn + 1, 3];
        }
    }

    public class CnvEvaluator
    {
        private readonly CNVChecker _cnvChecker;
        #region Members
        private const int MaxCn = 5; // Currently the max copynum is 5 

        #endregion

        public CnvEvaluator(CNVChecker cnvChecker)
        {
            _cnvChecker = cnvChecker;
        }

        public void ComputeAccuracy(Dictionary<string, List<CNInterval>> knownCN, string cnvCallsPath, string outputPath, bool includePassingOnly, EvaluateCnvOptions options, Dictionary<string, List<CnvCall>> calls)
        {
            // Make a note of how many bases in the truth set are not *actually* considered to be known bases, using
            // the "cnaqc" exclusion set:
            _cnvChecker.InitializeIntervalMetrics(knownCN);
            bool regionsOfInterest = _cnvChecker.RegionsOfInterest != null;
            var baseCounters = new List<BaseCounter> { new BaseCounter(MaxCn, 0, Int32.MaxValue, regionsOfInterest) };
            if (options.SplitBySize)
            {
                baseCounters.Add(new BaseCounter(MaxCn, 0, 4999, regionsOfInterest));
                baseCounters.Add(new BaseCounter(MaxCn, 5000, 9999, regionsOfInterest));
                baseCounters.Add(new BaseCounter(MaxCn, 10000, 99999, regionsOfInterest));
                baseCounters.Add(new BaseCounter(MaxCn, 100000, 499999, regionsOfInterest));
                baseCounters.Add(new BaseCounter(MaxCn, 500000, int.MaxValue, regionsOfInterest));
            }

            _cnvChecker.CountExcludedBasesInTruthSetIntervals(knownCN);

            foreach (var baseCounter in baseCounters)
            {
                var metrics = CalculateMetrics(knownCN, calls, baseCounter, options.SkipDiploid, includePassingOnly);

                string fileName = $"{options.BaseFileName}";
                if (options.DQscoreThreshold.HasValue)
                {
                    fileName += "_denovo";
                }
                if (baseCounter.MinSize != 0 || baseCounter.MaxSize != int.MaxValue)
                {
                    fileName += $"_{Math.Round(baseCounter.MinSize / 1000.0)}kb";
                    fileName += baseCounter.MaxSize == int.MaxValue ? "+" : $"_{ Math.Round(baseCounter.MaxSize / 1000.0)}kb";
                }
                fileName += ".txt";
                var outputDir = new DirectoryLocation(outputPath);
                outputDir.Create();
                var outputFile = outputDir.GetFileLocation(fileName);
                using (FileStream stream = new FileStream(outputFile.FullName, includePassingOnly ?
                FileMode.Create : FileMode.Append, FileAccess.Write))
                using (StreamWriter outputWriter = new StreamWriter(stream))
                {
                    outputWriter.NewLine = "\n";
                    WriteResults(cnvCallsPath, outputWriter, baseCounter, includePassingOnly, metrics);
                }
            }
        }

        public MetricsCalculator CalculateMetrics(Dictionary<string, List<CNInterval>> knownCN, Dictionary<string, List<CnvCall>> calls, 
            BaseCounter baseCounter, bool optionsSkipDiploid, bool includePassingOnly)
        {
            calls.Values.SelectMany(x => x).ForEach(call =>
            {
                if (!(call.IsAltVariant && call.Length >= baseCounter.MinSize && call.Length <= baseCounter.MaxSize))
                    return;
                if (includePassingOnly && !call.PassFilter)
                    return;
                baseCounter.TotalVariantBases += call.Length;
                baseCounter.TotalVariants++;
            });

            foreach (CNInterval interval in knownCN.Values.SelectMany(x=>x))
            {
                if (!(interval.Length >= baseCounter.MinSize && interval.Length <= baseCounter.MaxSize)) continue;
                int nonOverlapBases = interval.Length;
                int totalOverlapBases = 0;
                int excludeIntervalBases = 0;
                var totalIntervalRefPloidy = new List<(int ploidy, int length)>();
                string chromosome = interval.Chromosome;
                if (!calls.ContainsKey(chromosome)) chromosome = chromosome.Replace("chr", "");
                if (!calls.ContainsKey(chromosome)) chromosome = "chr" + chromosome;
                if (!calls.ContainsKey(chromosome))
                {
                    Console.WriteLine($"Error: Skipping truth variant for chromosome {interval.Chromosome} with no Canvas calls");
                    continue;
                }

                int knownCn = interval.Cn;
                if (knownCn > MaxCn) knownCn = MaxCn;

                foreach (CnvCall call in calls[chromosome])
                {
                    int CN = call.CN;
                    if (CN < 0 || call.End < 0) continue; // Not a CNV call, apparently
                    if (call.AltAllele == "." && optionsSkipDiploid) continue;

                    if (CN > MaxCn) CN = MaxCn;
                    string chr = call.Chr;
                    int overlapStart = Math.Max(call.Start, interval.Start);
                    int overlapEnd = Math.Min(call.End, interval.End);
                    if (overlapStart >= overlapEnd) continue;
                    int overlapBases = overlapEnd - overlapStart;
                    // We've got an overlap interval.  Kill off some bases from this interval, if it happens
                    // to overlap with an excluded interval:
                    if (_cnvChecker.ExcludeIntervals.ContainsKey(chr))
                    {
                        foreach (CNInterval excludeInterval in _cnvChecker.ExcludeIntervals[chr])
                        {
                            int excludeOverlapStart = Math.Max(excludeInterval.Start, overlapStart);
                            int excludeOverlapEnd = Math.Min(excludeInterval.End, overlapEnd);
                            if (excludeOverlapStart >= excludeOverlapEnd) continue;
                            excludeIntervalBases += excludeOverlapEnd - excludeOverlapStart;
                            overlapBases -= excludeOverlapEnd - excludeOverlapStart;
                            // if majority of the region is in exclude intervals, don't consider any overlap
                            if (overlapBases / Math.Max(excludeOverlapEnd - excludeOverlapStart, 1) < 0.1)
                                overlapBases = 0;
                        }
                    }

                    totalIntervalRefPloidy.Add((call.RefPloidy, overlapBases));

                    if (!call.PassFilter && includePassingOnly && knownCn != call.RefPloidy)
                        // assign no call (CN=ploidy) by default
                        baseCounter.BaseCount[knownCn, call.RefPloidy, call.RefPloidy] += overlapBases;
                    else
                    {
                        totalOverlapBases += overlapBases;
                        baseCounter.BaseCount[knownCn, CN, call.RefPloidy] += overlapBases;
                    }

                    interval.BasesCovered += overlapBases;

                    if (knownCn == CN)
                        interval.BasesCalledCorrectly += overlapBases;
                    else
                        interval.BasesCalledIncorrectly += overlapBases;

                    if (_cnvChecker.RegionsOfInterest == null ||
                        !_cnvChecker.RegionsOfInterest.ContainsKey(chr)) continue;

                    foreach (CNInterval roiInterval in _cnvChecker.RegionsOfInterest[chr])
                    {
                        int roiOverlapStart = Math.Max(roiInterval.Start, overlapStart);
                        int roiOverlapEnd = Math.Min(roiInterval.End, overlapEnd);
                        if (roiOverlapStart >= roiOverlapEnd) continue;
                        int roiOverlapBases = roiOverlapEnd - roiOverlapStart;
                        if (!call.PassFilter && includePassingOnly)
                            // assign no call (CN=ploidy) by default
                            baseCounter.RoiBaseCount[knownCn, call.RefPloidy, call.RefPloidy] += roiOverlapBases;
                        else
                            baseCounter.RoiBaseCount[knownCn, CN, call.RefPloidy] += roiOverlapBases;
                    }
                }

                nonOverlapBases -= (totalOverlapBases + excludeIntervalBases);
                if (totalIntervalRefPloidy.Empty())
                {
                    Console.WriteLine($"Error: Truth variant {interval.Chromosome}:{interval.Start}-{interval.End} with no overlapping " +
                                      $"Canvas calls. Ploidy cannot be determined!");
                }
                int ploidy = Convert.ToInt32(Math.Round(Utilities.WeightedMean(totalIntervalRefPloidy.Select(x => (double) x.ploidy).ToList(),
                    totalIntervalRefPloidy.Select(x => (double) Math.Max(x.length,1)).ToList())));
                interval.ReferenceCopyNumber = ploidy;
                if (nonOverlapBases < 0)
                {
                    throw new InvalidDataException($"Truth variant {interval.Chromosome}:{interval.Start}-{interval.End} has negative non-overlap bases");
                }
                baseCounter.BaseCount[knownCn, ploidy, ploidy] += nonOverlapBases;

            }

            CalculateMedianAndMeanAccuracies(baseCounter, knownCN);
            var allIntervals = knownCN.SelectMany(kvp => kvp.Value).ToList();

            // find truth interval with highest number of false negatives (hurts recall)
            var variantIntervals = allIntervals.Where(interval => interval.Cn != interval.ReferenceCopyNumber).ToList();
            if (variantIntervals.Any())
            {
                var intervalMaxFalseNegatives = variantIntervals.MaxBy(interval => interval.BasesNotCalled + interval.BasesCalledIncorrectly);
                Console.WriteLine($"Truth interval with most false negatives (hurts recall): {intervalMaxFalseNegatives}");
            }
            // find truth interval with highest number of false positive (hurts precision)
            var refIntervals = allIntervals.Where(interval => interval.Cn == interval.ReferenceCopyNumber).ToList();
            if (refIntervals.Any())
            {
                var intervalMaxFalsePositives = refIntervals.MaxBy(interval => interval.BasesCalledIncorrectly);
                Console.WriteLine($"Truth interval with most false positives (hurts precision): {intervalMaxFalsePositives}");
            }
            return MetricsCalculator.CalculateMetrics(baseCounter, MaxCn, 2);

        }

        /// <summary>
        /// For each CNV calls in the truth set, compute the fraction of bases assigned correct copy number
        /// </summary>
        /// <param name="baseCounter"></param>
        private void CalculateMedianAndMeanAccuracies(BaseCounter baseCounter, Dictionary<string, List<CNInterval>> knownCN)
        {
            baseCounter.MeanAccuracy = 0;
            baseCounter.MedianAccuracy = 0;
            var eventAccuracies = new List<double>();
            foreach (string chr in knownCN.Keys)
            {
                foreach (var interval in knownCN[chr])
                {
                    if (interval.Cn == interval.ReferenceCopyNumber) continue;
                    int basecount = interval.Length - interval.BasesExcluded;
                    if (basecount <= 0) continue;
                    double accuracy = interval.BasesCalledCorrectly / (double)basecount;
                    eventAccuracies.Add(accuracy);
                    baseCounter.MeanAccuracy += accuracy;
                    //Console.WriteLine("{0}\t{1:F4}", interval.End - interval.Start, accuracy);
                }
            }
            eventAccuracies.Sort();
            baseCounter.MeanAccuracy /= Math.Max(1, eventAccuracies.Count);
            baseCounter.MedianAccuracy = double.NaN;
            if (eventAccuracies.Count > 0)
                baseCounter.MedianAccuracy = eventAccuracies[eventAccuracies.Count / 2];
            Console.WriteLine($"Event-level accuracy mean {baseCounter.MeanAccuracy:F4} median {baseCounter.MedianAccuracy:F4}" +
                              $" for variants sizes {baseCounter.MinSize} to {baseCounter.MaxSize}");
        }

        private void WriteResults(string cnvCallsPath, StreamWriter outputWriter, BaseCounter baseCounter, bool includePassingOnly, MetricsCalculator metrics)
        {

            // load and append VCF header information 
            _cnvChecker.HandleVcfHeaderInfo(outputWriter, new FileLocation(cnvCallsPath));
            // Report stats:
            outputWriter.WriteLine(includePassingOnly ? "Results for PASSing variants" : "Results for all variants");
            outputWriter.WriteLine("Accuracy\t{0:F4}", metrics.Accuracy);
            // SK: I felt the direction based performance metrices make more sense
            outputWriter.WriteLine("DirectionAccuracy\t{0:F4}", metrics.DirectionAccuracy);
            outputWriter.WriteLine("F-score\t{0:F4}", metrics.F1Score);
            outputWriter.WriteLine("Recall\t{0:F4}", metrics.Recall);
            outputWriter.WriteLine("DirectionRecall\t{0:F4}", metrics.DirectionRecall);
            outputWriter.WriteLine("Precision\t{0:F4}", metrics.Precision);
            outputWriter.WriteLine("DirectionPrecision\t{0:F4}", metrics.DirectionPrecision);
            outputWriter.WriteLine("GainRecall\t{0:F4}", metrics.GainRecall);
            outputWriter.WriteLine("GainDirectionRecall\t{0:F4}", metrics.GainDirectionRecall);
            outputWriter.WriteLine("GainPrecision\t{0:F4}", metrics.GainPrecision);
            outputWriter.WriteLine("GainDirectionPrecision\t{0:F4}", metrics.GainDirectionPrecision);
            outputWriter.WriteLine("LossRecall\t{0:F4}", metrics.LossRecall);
            outputWriter.WriteLine("LossDirectionRecall\t{0:F4}", metrics.LossRecall);
            outputWriter.WriteLine("LossPrecision\t{0:F4}", metrics.LossPrecision);
            outputWriter.WriteLine("LossDirectionPrecision\t{0:F4}", metrics.LossDirectionPrecision);
            outputWriter.WriteLine("MeanEventAccuracy\t{0:F4}", 100 * baseCounter.MeanAccuracy);
            outputWriter.WriteLine("MedianEventAccuracy\t{0:F4}", 100 * baseCounter.MedianAccuracy);
            outputWriter.WriteLine("VariantEventsCalled\t{0}", baseCounter.TotalVariants);
            outputWriter.WriteLine("VariantBasesCalled\t{0}", baseCounter.TotalVariantBases);
            if (baseCounter.RoiBaseCount != null && metrics.RoiBases > 0)
            {
                outputWriter.WriteLine("ROIAccuracy\t{0:F4}", metrics.ROIAccuracy);
                outputWriter.WriteLine("ROIDirectionAccuracy\t{0:F4}", metrics.ROIDirectionAccuracy);
            }
            // to separate passing and all variant results
            outputWriter.WriteLine();
        }
    }
}