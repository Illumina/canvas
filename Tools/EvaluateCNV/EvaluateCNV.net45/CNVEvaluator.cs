using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using CanvasCommon;
using Illumina.Common;

namespace EvaluateCNV
{
    class CNVEvaluator
    {
        private CNVChecker _cnvChecker;
        #region Members
        static int maxCN = 10;
        readonly long[,] BaseCount = new long[maxCN + 1, maxCN + 1];
        readonly long[,] ROIBaseCount = new long[maxCN + 1, maxCN + 1];
        readonly long[,] BaseCountUnder5kb = new long[maxCN + 1, maxCN + 1];
        readonly long[,] BaseCount5kb10 = new long[maxCN + 1, maxCN + 1];
        readonly long[,] BaseCount10kb100kb = new long[maxCN + 1, maxCN + 1];
        readonly long[,] BaseCount100kb500 = new long[maxCN + 1, maxCN + 1];
        readonly long[,] BaseCountOver500kb = new long[maxCN + 1, maxCN + 1];
        int totalVariants = 0;
        long totalVariantBases = 0;
        #endregion

        public CNVEvaluator(CNVChecker cnvChecker)
        {
            _cnvChecker = cnvChecker;
        }

        public void ComputeAccuracy(string truthSetPath, string cnvCallsPath, string outputPath, PloidyInfo ploidyInfo, bool includePassingOnly, bool splitBySize)
        {
            // Make a note of how many bases in the truth set are not *actually* considered to be known bases, using
            // the "cnaqc" exclusion set:
            _cnvChecker.CountExcludedBasesInTruthSetIntervals();
            if (_cnvChecker.DQscoreThreshold.HasValue && !Path.GetFileName(cnvCallsPath).ToLower().Contains("vcf"))
                throw new ArgumentException("CNV.vcf must be in a vcf format when --dqscore option is used");
            IEnumerable<CNVCall> calls = Path.GetFileName(cnvCallsPath).ToLower().Contains("vcf")
                ? _cnvChecker.GetCnvCallsFromVcf(cnvCallsPath, includePassingOnly)
                : _cnvChecker.GetCnvCallsFromBed(cnvCallsPath);

            ploidyInfo.MakeChromsomeNameAgnosticWithAllChromosomes(calls.Select(call => call.Chr));
            foreach (CNVCall call in calls)
            {
                int CN = call.CN;
                if (CN < 0 || call.End < 0) continue; // Not a CNV call, apparently

                int basesOverlappingPloidyRegion = 0;
                int variantBasesOverlappingPloidyRegion = 0;
                foreach (PloidyInterval ploidyInterval in ploidyInfo.PloidyByChromosome[call.Chr])
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
                if (variantBasesOverlappingPloidyRegion > 0 ||
                    (CN != 2 && variantBasesOverlappingPloidyRegion < call.Length))
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
                    BaseCount[knownCN, CN] += overlapBases;
                    if (totalVariantBases < 5000)
                        BaseCountUnder5kb[knownCN, CN] += overlapBases;
                    else if (totalVariantBases >= 5000 && totalVariantBases <10000)
                        BaseCount5kb10[knownCN, CN] += overlapBases;
                    else if (totalVariantBases >= 10000 && totalVariantBases < 100000)
                        BaseCount10kb100kb[knownCN, CN] += overlapBases;
                    else if (totalVariantBases >= 100000 && totalVariantBases < 500000)
                        BaseCount100kb500[knownCN, CN] += overlapBases;
                    else if (totalVariantBases >= 500000)
                        BaseCountOver500kb[knownCN, CN] += overlapBases;

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
                            ROIBaseCount[knownCN, CN] += roiOverlapBases;
                        }
                    }
                }
            }

            // For each CNV calls in the truth set, compute the fraction of bases assigned correct copy number:
            List<double> eventAccuracies = new List<double>();
            double meanAccuracy = 0;
            foreach (string chr in _cnvChecker.KnownCN.Keys)
            {
                foreach (CNInterval interval in _cnvChecker.KnownCN[chr])
                {
                    if (interval.CN == 2) continue;
                    int baseCount = interval.Length - interval.BasesExcluded;
                    if (baseCount <= 0) continue;
                    double accuracy = interval.BasesCalledCorrectly / (double)baseCount;
                    eventAccuracies.Add(accuracy);
                    meanAccuracy += accuracy;
                    //Console.WriteLine("{0}\t{1:F4}", interval.End - interval.Start, accuracy);
                }
            }
            eventAccuracies.Sort();
            meanAccuracy /= Math.Max(1, eventAccuracies.Count);
            double medianAccuracy = double.NaN;
            if (eventAccuracies.Count > 0)
                medianAccuracy = eventAccuracies[eventAccuracies.Count / 2];
            Console.WriteLine("Event-level accuracy mean {0:F4} median {1:F4}", meanAccuracy, medianAccuracy);

            IEnumerable<CNInterval> allIntervals = _cnvChecker.KnownCN.SelectMany(kvp => kvp.Value);

            // find truth interval with highest number of false negatives (hurts recall)
            IEnumerable<CNInterval> variantIntervals = allIntervals.Where(interval => interval.CN != interval.ReferenceCopyNumber);
            CNInterval intervalMaxFalseNegatives = variantIntervals.MaxBy(interval => interval.BasesNotCalled + interval.BasesCalledIncorrectly);
            Console.WriteLine($"Truth interval with most false negatives (hurts recall): {intervalMaxFalseNegatives}");

            // find truth interval with highest number of false positive (hurts precision)
            IEnumerable<CNInterval> refIntervals = allIntervals.Where(interval => interval.CN == interval.ReferenceCopyNumber);
            CNInterval intervalMaxFalsePositives = refIntervals.MaxBy(interval => interval.BasesCalledIncorrectly);
            Console.WriteLine($"Truth interval with most false positives (hurts precision): {intervalMaxFalsePositives}");

            using (FileStream stream = new FileStream(Path.Combine(outputPath, "EvaluateCNV.txt"), FileMode.Create, FileAccess.Write))
            using (StreamWriter outputWriter = new StreamWriter(stream))
            {
                WriteResults(truthSetPath, cnvCallsPath, outputWriter, includePassingOnly, meanAccuracy, medianAccuracy, 
                    BaseCount, ROIBaseCount, totalVariants, totalVariantBases);
            }

            if (splitBySize)
            {
                using (FileStream stream = new FileStream(Path.Combine(outputPath, "EvaluateCNV_Under5kb.txt"), FileMode.Create, FileAccess.Write))
                using (StreamWriter outputWriter = new StreamWriter(stream))
                {
                    WriteResults(truthSetPath, cnvCallsPath, outputWriter, includePassingOnly, meanAccuracy, medianAccuracy,
                        BaseCountUnder5kb, null, totalVariants, totalVariantBases);
                }
                using (FileStream stream = new FileStream(Path.Combine(outputPath, "EvaluateCNV_5kb10.txt"), FileMode.Create, FileAccess.Write))
                using (StreamWriter outputWriter = new StreamWriter(stream))
                {
                    WriteResults(truthSetPath, cnvCallsPath, outputWriter, includePassingOnly, meanAccuracy, medianAccuracy,
                        BaseCount5kb10, null, totalVariants, totalVariantBases);
                }
                using (FileStream stream = new FileStream(Path.Combine(outputPath, "EvaluateCNV_10kb100kb.txt"), FileMode.Create, FileAccess.Write))
                using (StreamWriter outputWriter = new StreamWriter(stream))
                {
                    WriteResults(truthSetPath, cnvCallsPath, outputWriter, includePassingOnly, meanAccuracy, medianAccuracy,
                        BaseCount10kb100kb, null, totalVariants, totalVariantBases);
                }
                using (FileStream stream = new FileStream(Path.Combine(outputPath, "EvaluateCNV_100kb500.txt"), FileMode.Create, FileAccess.Write))
                using (StreamWriter outputWriter = new StreamWriter(stream))
                {
                    WriteResults(truthSetPath, cnvCallsPath, outputWriter, includePassingOnly, meanAccuracy, medianAccuracy,
                        BaseCount100kb500, null, totalVariants, totalVariantBases);
                }
            }
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
                        if ((trueCN < 2 && callCN < 2) || (trueCN == 2 && callCN == 2) || (trueCN > 2 && callCN > 2))
                            ROIBasesCorrectDirection += bases;
                    }
                }
            }

            // Report stats:
            if (includePassingOnly)
            {
                outputWriter.WriteLine("Results for PASSing variants");
            }
            else
            {
                outputWriter.WriteLine("Results for all variants");
            }
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