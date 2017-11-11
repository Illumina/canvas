using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using CanvasCommon;
using Illumina.Common;
using Illumina.Common.FileSystem;

namespace EvaluateCNV
{
    internal class BaseCounter
    {
        public int TotalVariants { get; set; }
        public int TotalVariantBases { get; set; }
        public double MeanAccuracy { get; set; }
        public double MedianAccuracy { get; set; }
        public int MinSize { get; }
        public int MaxSize { get; }
        public long[,] BaseCount;
        public long[,] RoiBaseCount;


        public BaseCounter(int maxCn, int minSize, int maxSize, bool hasRoi = false)
        {
            MinSize = minSize;
            MaxSize = maxSize;
            BaseCount = new long[maxCn + 1, maxCn + 1];
            if (hasRoi)
                RoiBaseCount = new long[maxCn + 1, maxCn + 1];
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
                if (!(call.Length >= baseCounter.MinSize && call.Length <= baseCounter.MaxSize)) continue;

                int basesOverlappingPloidyRegion = 0;
                int variantBasesOverlappingPloidyRegion = 0;
                foreach (var ploidyInterval in ploidyInfo.PloidyByChromosome[call.Chr])
                {
                    int overlap = call.Overlap(ploidyInterval);
                    basesOverlappingPloidyRegion += overlap;
                    if (CN != ploidyInterval.Ploidy)
                        variantBasesOverlappingPloidyRegion += overlap;
                }
                baseCounter.TotalVariantBases += variantBasesOverlappingPloidyRegion;
                if (CN != 2)
                {
                    baseCounter.TotalVariantBases += call.Length - basesOverlappingPloidyRegion;
                }
                if (variantBasesOverlappingPloidyRegion > 0 || CN != 2 && variantBasesOverlappingPloidyRegion < call.Length)
                {
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
                    baseCounter.BaseCount[knownCn, CN] += overlapBases;
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
                            baseCounter.RoiBaseCount[knownCn, CN] += roiOverlapBases;
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
            // Compute overall stats:
            long totalBases = 0;
            long totalBasesRight = 0;
            long totalBasesRightDirection = 0;

            long isGainBases = 0;
            long callGainBases = 0;
            long isGainBasesCorrect = 0;
            long isGainBasesCorrectDirection = 0;

            long isLossBases = 0;
            long callLossBases = 0;
            long isLossBasesCorrect = 0;
            long isLossBasesCorrectDirection = 0;
            for (int trueCn = 0; trueCn <= MaxCn; trueCn++)
            {
                for (int callCn = 0; callCn <= MaxCn; callCn++)
                {
                    long bases = baseCounter.BaseCount[trueCn, callCn];
                    totalBases += bases;
                    if (trueCn == callCn) totalBasesRight += bases;
                    if (trueCn < 2 && callCn < 2 || trueCn == 2 && callCn == 2 || trueCn > 2 && callCn > 2)
                        totalBasesRightDirection += bases;
                    if (trueCn < 2) isLossBases += bases;
                    if (trueCn > 2) isGainBases += bases;
                    if (callCn < 2) callLossBases += bases;
                    if (callCn > 2) callGainBases += bases;
                    if (trueCn == callCn && trueCn < 2) isLossBasesCorrect += bases;
                    if (trueCn == callCn && trueCn > 2) isGainBasesCorrect += bases;
                    if (trueCn > 2 && callCn > 2) isGainBasesCorrectDirection += bases;
                    if (trueCn < 2 && callCn < 2) isLossBasesCorrectDirection += bases;
                }
            }

            // Compute ROI stats:
            long roiBases = 0;
            long roiBasesCorrect = 0;
            long roiBasesCorrectDirection = 0;
            if (baseCounter.RoiBaseCount != null)
            {
                for (int trueCn = 0; trueCn <= MaxCn; trueCn++)
                {
                    for (int callCn = 0; callCn <= MaxCn; callCn++)
                    {
                        long bases = baseCounter.RoiBaseCount[trueCn, callCn];
                        roiBases += bases;
                        if (trueCn == callCn) roiBasesCorrect += bases;
                        if (trueCn < 2 && callCn < 2 || trueCn == 2 && callCn == 2 || trueCn > 2 && callCn > 2)
                            roiBasesCorrectDirection += bases;
                    }
                }
            }

            // load and append VCF header information 
            _cnvChecker.HandleVcfHeaderInfo(outputWriter, new FileLocation(cnvCallsPath));

            // Report stats:
            var precision = (isGainBasesCorrect + isLossBasesCorrect) / (double)(callGainBases + callLossBases);
            var recall = (isGainBasesCorrect + isLossBasesCorrect) / (double)(isGainBases + isLossBases);
            var f1Score = 2 * precision * recall / (precision + recall);
            outputWriter.WriteLine(includePassingOnly ? "Results for PASSing variants" : "Results for all variants");
            outputWriter.WriteLine("Accuracy\t{0:F4}", 100 * totalBasesRight / (double)totalBases);
            // SK: I felt the direction based performance metrices make more sense
            outputWriter.WriteLine("DirectionAccuracy\t{0:F4}", 100 * totalBasesRightDirection / (double)totalBases);
            outputWriter.WriteLine("F-score\t{0:F4}", f1Score);
            outputWriter.WriteLine("Recall\t{0:F4}", 100 * recall);
            outputWriter.WriteLine("DirectionRecall\t{0:F4}", 100 * (isGainBasesCorrectDirection + isLossBasesCorrectDirection) / (double)(isGainBases + isLossBases));
            outputWriter.WriteLine("Precision\t{0:F4}", 100 * precision);
            outputWriter.WriteLine("DirectionPrecision\t{0:F4}",
                100 * (isGainBasesCorrectDirection + isLossBasesCorrectDirection) / (double)(callGainBases + callLossBases));
            outputWriter.WriteLine("GainRecall\t{0:F4}", 100 * (isGainBasesCorrect) / (double)(isGainBases));
            outputWriter.WriteLine("GainDirectionRecall\t{0:F4}", 100 * (isGainBasesCorrectDirection) / (double)(isGainBases));
            outputWriter.WriteLine("GainPrecision\t{0:F4}", 100 * (isGainBasesCorrect) / (double)(callGainBases));
            outputWriter.WriteLine("GainDirectionPrecision\t{0:F4}",
                100 * (isGainBasesCorrectDirection) / (double)(callGainBases));
            outputWriter.WriteLine("LossRecall\t{0:F4}", 100 * (isLossBasesCorrect) / (double)(isLossBases));
            outputWriter.WriteLine("LossDirectionRecall\t{0:F4}", 100 * (isLossBasesCorrectDirection) / (double)(isLossBases));
            outputWriter.WriteLine("LossPrecision\t{0:F4}", 100 * (isLossBasesCorrect) / (double)(callLossBases));
            outputWriter.WriteLine("LossDirectionPrecision\t{0:F4}",
                100 * (isLossBasesCorrectDirection) / (double)(callLossBases));
            outputWriter.WriteLine("MeanEventAccuracy\t{0:F4}", 100 * baseCounter.MeanAccuracy);
            outputWriter.WriteLine("MedianEventAccuracy\t{0:F4}", 100 * baseCounter.MedianAccuracy);
            outputWriter.WriteLine("VariantEventsCalled\t{0}", baseCounter.TotalVariants);
            outputWriter.WriteLine("VariantBasesCalled\t{0}", baseCounter.TotalVariantBases);
            if (baseCounter.RoiBaseCount != null && roiBases > 0)
            {
                outputWriter.WriteLine("ROIAccuracy\t{0:F4}", 100 * roiBasesCorrect / (double)roiBases);
                outputWriter.WriteLine("ROIDirectionAccuracy\t{0:F4}", 100 * roiBasesCorrectDirection / (double)roiBases);
            }
            // to separate passing and all variant results
            outputWriter.WriteLine();
        }
    }
}