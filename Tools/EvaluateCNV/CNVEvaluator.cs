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
    internal class BaseCounter
    {
        public int TotalVariants { get; set; }
        public long TotalVariantBases { get; set; }
        public double MeanAccuracy { get; set; }
        public double MedianAccuracy { get; set; }
        public int MinSize { get; }
        public int MaxSize { get; }
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

    class CnvEvaluator
    {
        private readonly CNVChecker _cnvChecker;
        #region Members
        private const int MaxCn = 10; // Currently the max copynum is 5 

        #endregion

        public CnvEvaluator(CNVChecker cnvChecker)
        {
            _cnvChecker = cnvChecker;
        }

        public void ComputeAccuracy(string truthSetPath, string cnvCallsPath, string outputPath, PloidyInfo ploidyInfo, bool includePassingOnly,
            EvaluateCnvOptions options)
        {
            // Make a note of how many bases in the truth set are not *actually* considered to be known bases, using
            // the "cnaqc" exclusion set:
            _cnvChecker.InitializeIntervalMetrics();
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

            _cnvChecker.CountExcludedBasesInTruthSetIntervals();
            if (_cnvChecker.DQscoreThreshold.HasValue && !Path.GetFileName(cnvCallsPath).ToLower().Contains("vcf"))
                throw new ArgumentException("CNV.vcf must be in a vcf format when --dqscore option is used");
            var calls = _cnvChecker.GetCnvCallsFromVcf(cnvCallsPath, includePassingOnly);

            foreach (var baseCounter in baseCounters)
            {
                CalculateMetrics(ploidyInfo, calls, baseCounter, options.SkipDiploid);
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
                    WriteResults(cnvCallsPath, outputWriter, baseCounter, includePassingOnly);
                }
            }
        }

        private void CalculateMetrics(PloidyInfo ploidyInfo, IEnumerable<CnvCall> calls, BaseCounter baseCounter, bool optionsSkipDiploid)
        {
            ploidyInfo.MakeChromsomeNameAgnosticWithAllChromosomes(calls.Select(call => call.Chr));
            foreach (CnvCall call in calls)
            {
                int CN = call.CN;
                if (CN < 0 || call.End < 0) continue; // Not a CNV call, apparently
                if (call.AltAllele == "." && optionsSkipDiploid) continue;

                if (call.IsAltVariant)
                {
                    baseCounter.TotalVariantBases += call.Length;
                    baseCounter.TotalVariants++;
                }

                if (CN > MaxCn) CN = MaxCn;
                string chr = call.Chr;
                if (!_cnvChecker.KnownCn.ContainsKey(chr)) chr = call.Chr.Replace("chr", "");
                if (!_cnvChecker.KnownCn.ContainsKey(chr)) chr = "chr" + call.Chr;
                if (!_cnvChecker.KnownCn.ContainsKey(chr))
                {
                    Console.WriteLine("Error: Skipping variant call for chromosome {0} with no truth data", call.Chr);
                    continue;
                }
                foreach (CNInterval interval in _cnvChecker.KnownCn[chr])
                {
                    if (!(interval.Length >= baseCounter.MinSize && interval.Length <= baseCounter.MaxSize)) continue;

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
                            overlapBases -= (excludeOverlapEnd - excludeOverlapStart);
                        }
                    }

                    int knownCn = interval.Cn;
                    if (knownCn > MaxCn) knownCn = MaxCn;
                    // BasesCovered does not model ploidy-specific information, adjust CN ploidy instead 

                    baseCounter.BaseCount[knownCn, CN, call.RefPloidy] += overlapBases;
                    interval.BasesCovered += overlapBases;
                    if (knownCn == CN)
                    {
                        interval.BasesCalledCorrectly += overlapBases;
                    }
                    else
                    {
                        interval.BasesCalledIncorrectly += overlapBases;
                    }

                    if (_cnvChecker.RegionsOfInterest != null && _cnvChecker.RegionsOfInterest.ContainsKey(chr))
                    {
                        foreach (CNInterval roiInterval in _cnvChecker.RegionsOfInterest[chr])
                        {
                            int roiOverlapStart = Math.Max(roiInterval.Start, overlapStart);
                            int roiOverlapEnd = Math.Min(roiInterval.End, overlapEnd);
                            if (roiOverlapStart >= roiOverlapEnd) continue;
                            int roiOverlapBases = roiOverlapEnd - roiOverlapStart;
                            baseCounter.RoiBaseCount[knownCn, CN, call.RefPloidy] += roiOverlapBases;
                        }
                    }
                }
            }

            CalculateMedianAndMeanAccuracies(baseCounter);

            var allIntervals = _cnvChecker.KnownCn.SelectMany(kvp => kvp.Value).ToList();

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
        }

        /// <summary>
        /// For each CNV calls in the truth set, compute the fraction of bases assigned correct copy number
        /// </summary>
        /// <param name="baseCounter"></param>
        private void CalculateMedianAndMeanAccuracies(BaseCounter baseCounter)
        {
            baseCounter.MeanAccuracy = 0;
            baseCounter.MedianAccuracy = 0;
            var eventAccuracies = new List<double>();
            foreach (string chr in _cnvChecker.KnownCn.Keys)
            {
                foreach (var interval in _cnvChecker.KnownCn[chr])
                {
                    if (interval.Cn == 2) continue;
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

        private void WriteResults(string cnvCallsPath, StreamWriter outputWriter, BaseCounter baseCounter, bool includePassingOnly)
        {
           
            // load and append VCF header information 
            _cnvChecker.HandleVcfHeaderInfo(outputWriter, new FileLocation(cnvCallsPath));
            var metrics = MetricsCalculator.CalculateMetrics(baseCounter, MaxCn, 2);

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