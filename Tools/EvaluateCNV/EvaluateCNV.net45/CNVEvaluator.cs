using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using CanvasCommon;
using Illumina.Common;

namespace EvaluateCNV
{
    internal class BaseCounter
    {
        public int? MinSize { get; }
        public int? MaxSize { get; }
        public bool HasRoi { get; }
        public long[,] BaseCount;

        public BaseCounter(int maxCN, int? minSize = null, int? maxSize = null, bool hasRoi = false)
        {
            MinSize = minSize;
            MaxSize = maxSize;
            HasRoi = hasRoi;
            BaseCount = new long[maxCN + 1, maxCN + 1];
        }
    }

    class CNVEvaluator
    {
        private CNVChecker _cnvChecker;
        #region Members
        static int maxCN = 10;
        #endregion

        public CNVEvaluator(CNVChecker cnvChecker)
        {
            _cnvChecker = cnvChecker;
        }

        public void ComputeAccuracy(string truthSetPath, string cnvCallsPath, string outputPath, PloidyInfo ploidyInfo, bool includePassingOnly,
            EvaluateCnvOptions options)
        {
            // Make a note of how many bases in the truth set are not *actually* considered to be known bases, using
            // the "cnaqc" exclusion set:


            var globalBaseCounter = new BaseCounter(maxCN);
            var roiBaseCounter = new BaseCounter(maxCN, minSize: null, maxSize: null, hasRoi: true);
            var sizeAwareBaseCounters = new List<BaseCounter>
            {
                new BaseCounter(maxCN, 0, 4999),
                new BaseCounter(maxCN, 5000, 9999),
                new BaseCounter(maxCN, 10000, 99999),
                new BaseCounter(maxCN, 100000, 499999),
                new BaseCounter(maxCN, 500000, int.MaxValue)
            };

            int totalVariants = 0;
            long totalVariantBases = 0;

            _cnvChecker.CountExcludedBasesInTruthSetIntervals();
            if (_cnvChecker.DQscoreThreshold.HasValue && !includePassingOnly)
                throw new ArgumentException("CNV.vcf must be in a vcf format when --dqscore option is used");
            var calls = _cnvChecker.GetCnvCallsFromVcf(cnvCallsPath, includePassingOnly);

            double medianAccuracy;
            var meanAccuracy = CalculateMetrics(ploidyInfo, calls, globalBaseCounter, sizeAwareBaseCounters, roiBaseCounter,
                options.SkipDiploid, ref totalVariantBases, ref totalVariants, out medianAccuracy);

            using (var stream = new FileStream(Path.Combine(outputPath, $"{options.BaseFileName}.txt"), includePassingOnly? 
                FileMode.Create : FileMode.Append, FileAccess.Write))
            using (StreamWriter outputWriter = new StreamWriter(stream))
            {
                WriteResults(truthSetPath, cnvCallsPath, outputWriter, includePassingOnly, meanAccuracy, medianAccuracy, 
                    globalBaseCounter.BaseCount, roiBaseCounter.BaseCount, totalVariants, totalVariantBases);
            }

       
            if (options.SplitBySize)
            {
                foreach (var baseCounter in sizeAwareBaseCounters)
                {
                    var fileName = $"{options.BaseFileName}_{Math.Round(baseCounter.MinSize.Value / 1000.0)}kb_" +
                                   $"{Math.Round(baseCounter.MaxSize.Value / 1000.0)}kb.txt";

                    using (FileStream stream = new FileStream(Path.Combine(outputPath, fileName), includePassingOnly ?
                        FileMode.Create : FileMode.Append, FileAccess.Write))
                    using (StreamWriter outputWriter = new StreamWriter(stream))
                    {
                        WriteResults(truthSetPath, cnvCallsPath, outputWriter, includePassingOnly, meanAccuracy, medianAccuracy,
                            baseCounter.BaseCount, null, totalVariants, totalVariantBases);
                    }
                }

            }
        }
            
        private double CalculateMetrics(PloidyInfo ploidyInfo, IEnumerable<CNVCall> calls, BaseCounter globalCounter, 
            List<BaseCounter> sizeAwareCounter, BaseCounter roiCounter, bool optionsSkipDiploid, ref long totalVariantBases, 
            ref int totalVariants, out double medianAccuracy)
        {
            ploidyInfo.MakeChromsomeNameAgnosticWithAllChromosomes(calls.Select(call => call.Chr));
            foreach (CNVCall call in calls)
            {
                int CN = call.CN;
                if (CN < 0 || call.End < 0) continue; // Not a CNV call, apparently
                if (call.AltAllele == "." && optionsSkipDiploid) continue;

                int basesOverlappingPloidyRegion = 0;
                int variantBasesOverlappingPloidyRegion = 0;
                foreach (var ploidyInterval in ploidyInfo.PloidyByChromosome[call.Chr])
                {
                    int overlap = call.Overlap(ploidyInterval);
                    basesOverlappingPloidyRegion += overlap;
                    if (CN != ploidyInterval.Ploidy)
                        variantBasesOverlappingPloidyRegion += overlap;
                }
                totalVariantBases += variantBasesOverlappingPloidyRegion;
                if (CN != 2)
                {
                    totalVariantBases += call.Length - basesOverlappingPloidyRegion;
                }
                if (variantBasesOverlappingPloidyRegion > 0 || CN != 2 && variantBasesOverlappingPloidyRegion < call.Length)
                {
                    totalVariants++;
                }

                if (CN > maxCN) CN = maxCN;
                string chr = call.Chr;
                if (!_cnvChecker.KnownCN.ContainsKey(chr)) chr = call.Chr.Replace("chr", "");
                if (!_cnvChecker.KnownCN.ContainsKey(chr)) chr = "chr" + call.Chr;
                if (!_cnvChecker.KnownCN.ContainsKey(chr))
                {
                    Console.WriteLine("Error: Skipping variant call for chromosome {0} with no truth data", call.Chr);
                    continue;
                }
                foreach (CNInterval interval in _cnvChecker.KnownCN[chr])
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

                    int knownCN = interval.CN;
                    if (knownCN > maxCN) knownCN = maxCN;
                    globalCounter.BaseCount[knownCN, CN] += overlapBases;
                    foreach (var baseCounter in sizeAwareCounter)
                    {
                        if (call.Length >= baseCounter.MinSize.Value && call.Length <= baseCounter.MaxSize.Value)
                            baseCounter.BaseCount[knownCN, CN] += overlapBases;
                    }

                    interval.BasesCovered += overlapBases;
                    if (knownCN == CN)
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
                            roiCounter.BaseCount[knownCN, CN] += roiOverlapBases;
                        }
                    }
                }
            }

            // For each CNV calls in the truth set, compute the fraction of bases assigned correct copy number:
            List<double> eventAccuracies = new List<double>();
            double meanAccuracy = 0;
            foreach (string chr in _cnvChecker.KnownCN.Keys)
            {
                foreach (var interval in _cnvChecker.KnownCN[chr])
                {
                    if (interval.CN == 2) continue;
                    int basecount = interval.Length - interval.BasesExcluded;
                    if (basecount <= 0) continue;
                    double accuracy = interval.BasesCalledCorrectly / (double) basecount;
                    eventAccuracies.Add(accuracy);
                    meanAccuracy += accuracy;
                    //Console.WriteLine("{0}\t{1:F4}", interval.End - interval.Start, accuracy);
                }
            }
            eventAccuracies.Sort();
            meanAccuracy /= Math.Max(1, eventAccuracies.Count);
            medianAccuracy = double.NaN;
            if (eventAccuracies.Count > 0)
                medianAccuracy = eventAccuracies[eventAccuracies.Count / 2];
            Console.WriteLine("Event-level accuracy mean {0:F4} median {1:F4}", meanAccuracy, medianAccuracy);

            var allIntervals = _cnvChecker.KnownCN.SelectMany(kvp => kvp.Value);

            // find truth interval with highest number of false negatives (hurts recall)
            var variantIntervals = allIntervals.Where(interval => interval.CN != interval.ReferenceCopyNumber);
            var intervalMaxFalseNegatives = variantIntervals.MaxBy(interval => interval.BasesNotCalled + interval.BasesCalledIncorrectly);
            Console.WriteLine($"Truth interval with most false negatives (hurts recall): {intervalMaxFalseNegatives}");

            // find truth interval with highest number of false positive (hurts precision)
            var refIntervals = allIntervals.Where(interval => interval.CN == interval.ReferenceCopyNumber);
            var intervalMaxFalsePositives = refIntervals.MaxBy(interval => interval.BasesCalledIncorrectly);
            Console.WriteLine($"Truth interval with most false positives (hurts precision): {intervalMaxFalsePositives}");
            return meanAccuracy;
        }

        static void WriteResults(string truthSetPath, string cnvCallsPath, StreamWriter outputWriter, bool includePassingOnly,
            double meanAccuracy, double medianAccuracy, long[,] baseCount, long[,] roiBaseCount, int totalVariants, long totalVariantBases)
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
            for (int trueCN = 0; trueCN <= maxCN; trueCN++)
            {
                for (int callCN = 0; callCN <= maxCN; callCN++)
                {
                    long bases = baseCount[trueCN, callCN];
                    totalBases += bases;
                    if (trueCN == callCN) totalBasesRight += bases;
                    if (trueCN < 2 && callCN < 2 || trueCN == 2 && callCN == 2 || trueCN > 2 && callCN > 2)
                        totalBasesRightDirection += bases;
                    if (trueCN < 2) isLossBases += bases;
                    if (trueCN > 2) isGainBases += bases;
                    if (callCN < 2) callLossBases += bases;
                    if (callCN > 2) callGainBases += bases;
                    if (trueCN == callCN && trueCN < 2) isLossBasesCorrect += bases;
                    if (trueCN == callCN && trueCN > 2) isGainBasesCorrect += bases;
                    if (trueCN > 2 && callCN > 2) isGainBasesCorrectDirection += bases;
                    if (trueCN < 2 && callCN < 2) isLossBasesCorrectDirection += bases;
                }
            }

            // Compute ROI stats:
            long ROIBases = 0;
            long ROIBasesCorrect = 0;
            long ROIBasesCorrectDirection = 0;
            if (roiBaseCount != null)
            {
                for (int trueCN = 0; trueCN <= maxCN; trueCN++)
                {
                    for (int callCN = 0; callCN <= maxCN; callCN++)
                    {
                        long bases = roiBaseCount[trueCN, callCN];
                        ROIBases += bases;
                        if (trueCN == callCN) ROIBasesCorrect += bases;
                        if (trueCN < 2 && callCN < 2 || trueCN == 2 && callCN == 2 || trueCN > 2 && callCN > 2)
                            ROIBasesCorrectDirection += bases;
                    }
                }
            }

            // Report stats:
            outputWriter.WriteLine(includePassingOnly ? "Results for PASSing variants" : "Results for all variants");

            outputWriter.WriteLine("TruthSet\t{0}", truthSetPath);
            outputWriter.WriteLine("CNVCalls\t{0}", cnvCallsPath);
            outputWriter.WriteLine("Accuracy\t{0:F4}", 100 * totalBasesRight / (double) totalBases);
            outputWriter.WriteLine("DirectionAccuracy\t{0:F4}", 100 * totalBasesRightDirection / (double) totalBases);
            outputWriter.WriteLine("Recall\t{0:F4}",
                100 * (isGainBasesCorrect + isLossBasesCorrect) / (double) (isGainBases + isLossBases));
            outputWriter.WriteLine("DirectionRecall\t{0:F4}",
                100 * (isGainBasesCorrectDirection + isLossBasesCorrectDirection) / (double) (isGainBases + isLossBases));
            outputWriter.WriteLine("Precision\t{0:F4}",
                100 * (isGainBasesCorrect + isLossBasesCorrect) / (double) (callGainBases + callLossBases));
            outputWriter.WriteLine("DirectionPrecision\t{0:F4}",
                100 * (isGainBasesCorrectDirection + isLossBasesCorrectDirection) / (double) (callGainBases + callLossBases));
            outputWriter.WriteLine("GainRecall\t{0:F4}", 100 * (isGainBasesCorrect) / (double) (isGainBases));
            outputWriter.WriteLine("GainDirectionRecall\t{0:F4}", 100 * (isGainBasesCorrectDirection) / (double) (isGainBases));
            outputWriter.WriteLine("GainPrecision\t{0:F4}", 100 * (isGainBasesCorrect) / (double) (callGainBases));
            outputWriter.WriteLine("GainDirectionPrecision\t{0:F4}",
                100 * (isGainBasesCorrectDirection) / (double) (callGainBases));
            outputWriter.WriteLine("LossRecall\t{0:F4}", 100 * (isLossBasesCorrect) / (double) (isLossBases));
            outputWriter.WriteLine("LossDirectionRecall\t{0:F4}", 100 * (isLossBasesCorrectDirection) / (double) (isLossBases));
            outputWriter.WriteLine("LossPrecision\t{0:F4}", 100 * (isLossBasesCorrect) / (double) (callLossBases));
            outputWriter.WriteLine("LossDirectionPrecision\t{0:F4}",
                100 * (isLossBasesCorrectDirection) / (double) (callLossBases));
            outputWriter.WriteLine("MeanEventAccuracy\t{0:F4}", 100 * meanAccuracy);
            outputWriter.WriteLine("MedianEventAccuracy\t{0:F4}", 100 * medianAccuracy);
            outputWriter.WriteLine("VariantEventsCalled\t{0}", totalVariants);
            outputWriter.WriteLine("VariantBasesCalled\t{0}", totalVariantBases);
            if (roiBaseCount != null && ROIBases > 0)
            {
                outputWriter.WriteLine("ROIAccuracy\t{0:F4}", 100 * ROIBasesCorrect / (double) ROIBases);
                outputWriter.WriteLine("ROIDirectionAccuracy\t{0:F4}", 100 * ROIBasesCorrectDirection / (double) ROIBases);
            }
            // to separate passing and all variant results
            outputWriter.WriteLine();
        }
    }
}