using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using CanvasCommon;
using Illumina.Common;
using Illumina.Common.FileSystem;
using Isas.SequencingFiles;

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
        public long[,] NoCalls; // 2D array stores knownCN and REF ploidy for bases where there are no calls


        public BaseCounter(int maxCn, int minSize, int maxSize, bool hasRoi = false)
        {
            MinSize = minSize;
            MaxSize = maxSize;
            BaseCount = new long[maxCn + 1, maxCn + 1, 3];
            NoCalls = new long[maxCn + 1, 3];
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
            bool regionsOfInterest = !_cnvChecker.RegionsOfInterest.Empty();
            var baseCounters = new List<BaseCounter> { new BaseCounter(MaxCn, 0, Int32.MaxValue, regionsOfInterest) };
            if (options.SplitBySize)
            {
                baseCounters.Add(new BaseCounter(MaxCn, 0, 4999, regionsOfInterest));
                baseCounters.Add(new BaseCounter(MaxCn, 5000, 9999, regionsOfInterest));
                baseCounters.Add(new BaseCounter(MaxCn, 10000, 99999, regionsOfInterest));
                baseCounters.Add(new BaseCounter(MaxCn, 100000, 499999, regionsOfInterest));
                baseCounters.Add(new BaseCounter(MaxCn, 500000, int.MaxValue, regionsOfInterest));
            }

            // not parallel here as parallelism will be attained at the level of regression workflow 
            _cnvChecker.CountExcludedBasesInTruthSetIntervals(knownCN);
            Dictionary<string, BitArray> referenceBases = null;
            if (options.KmerFa != null)
            {
                referenceBases = new Dictionary<string, BitArray>();
                foreach (var chr in knownCN.Keys)
                {
                    string chromReferenceBases = FastaLoader.LoadFastaSequence(options.KmerFa, chr);
                    var bitArrayBases = new BitArray(chromReferenceBases.Length);
                    // Mark which k-mers in the fasta file are unique. These are indicated by upper-case letters.
                    for (var i = 0; i < chromReferenceBases.Length; i++)
                    {
                        if (char.IsUpper(chromReferenceBases[i]))
                            bitArrayBases[i] = true;
                    }
                    referenceBases[chr] = bitArrayBases;
                }
            }

            foreach (var baseCounter in baseCounters)
            {
                _cnvChecker.InitializeIntervalMetrics(knownCN);
                var metrics = CalculateMetrics(knownCN, calls, baseCounter, options.SkipDiploid, includePassingOnly, referenceBases);

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
            BaseCounter baseCounter, bool optionsSkipDiploid, bool includePassingOnly, Dictionary<string, BitArray> kmerfa = null)
        {
            // string referenceBases = string.Empty;
            calls.Values.SelectMany(x => x).ForEach(call =>
            {
                if (!(call.IsAltVariant && call.Length >= baseCounter.MinSize && call.Length <= baseCounter.MaxSize))
                    return;
                if (includePassingOnly && !call.PassFilter)
                    return;
                baseCounter.TotalVariantBases += call.Length;
                baseCounter.TotalVariants++;
            });

            // skip truth interval that have >= 80% of unmapable bases
            // code is not parallel as this will be done by Regression workflow
            const double fractionUnmappableBases = 0.8;
            var filteredknownCn = new Dictionary<string, List<CNInterval>>();
            if (kmerfa != null)
            {
                foreach (var chromosome in knownCN.Keys)
                {
                    filteredknownCn[chromosome] = new List<CNInterval>();
                    foreach (var interval in knownCN[chromosome])
                    {
                        // always include REF intervals even if they are in unmappable regions
                        if (interval.Cn == interval.ReferenceCopyNumber)
                        {
                            filteredknownCn[chromosome].Add(interval);
                            continue;
                        }
                        var flaggedBasesCounter = 0;
                        for (var bp = interval.Start; bp < interval.End; bp++)
                        {
                            if (!kmerfa[chromosome][bp])
                            {
                                flaggedBasesCounter++;
                            }
                        }

                        if (flaggedBasesCounter / (double)interval.Length < fractionUnmappableBases)
                        {
                            filteredknownCn[chromosome].Add(interval);
                        }
                        else
                        {
                            Console.Error.WriteLine($"skipping truth interval {interval} with >= {fractionUnmappableBases} fraction of unmappable positions");
                        }
                    }
                }

            }
            else
            {
                filteredknownCn = knownCN;
            }

            foreach (CNInterval interval in filteredknownCn.Values.SelectMany(x => x))
            {
                if (!(interval.Length >= baseCounter.MinSize && interval.Length <= baseCounter.MaxSize)) continue;
                int nonOverlapBases = interval.Length;
                int nonOverlapRoiBases = 0;
                if (!_cnvChecker.RegionsOfInterest.Empty() &&
                    _cnvChecker.RegionsOfInterest.ContainsKey(interval.Chromosome))
                {

                    foreach (CNInterval roiInterval in _cnvChecker.RegionsOfInterest[interval.Chromosome])
                    {
                        int roiOverlapStart = Math.Max(roiInterval.Start, interval.Start);
                        int roiOverlapEnd = Math.Min(roiInterval.End, interval.End);
                        if (roiOverlapStart >= roiOverlapEnd) continue;
                        int roiOverlapBases = roiOverlapEnd - roiOverlapStart;
                        nonOverlapRoiBases -= roiOverlapBases;
                    }
                }
                int totalOverlapBases = 0;
                int totalRoiOverlapBases = 0;
                int excludeIntervalBases = 0;
                var totalIntervalRefPloidy = new List<(int ploidy, int length)>();
                string chromosome = interval.Chromosome;


                if (!calls.ContainsKey(chromosome)) chromosome = chromosome.Replace("chr", "");
                if (!calls.ContainsKey(chromosome)) chromosome = "chr" + chromosome;


                IEnumerable<CnvCall> callsThisChromosome;
                if (calls.ContainsKey(chromosome))
                {
                    callsThisChromosome = calls[chromosome];
                }
                else
                {
                    Console.Error.WriteLine($"Error: no Canvas calls for chromosome {interval.Chromosome} in truth file");
                    callsThisChromosome = Enumerable.Empty<CnvCall>();
                }
                int knownCn = interval.Cn;
                if (knownCn > MaxCn) knownCn = MaxCn;

                int thisIntervalBasesTruePositive = 0;
                int thisIntervalBasesTrueNegative = 0;
                int thisIntervalBasesFalsePositive = 0;
                int thisIntervalBasesFalseNegative = 0;
                int thisIntervalBasesNoCall = interval.Length;
                int thisIntervalBasesExcluded = 0;

                foreach (CnvCall call in callsThisChromosome)
                {
                    if (!call.RefPloidy.HasValue)
                        throw new IlluminaException($"Could not determine reference ploidy for call '{call}'. Please provide ploidy information via command line option.");
                    int refPloidy = interval.ReferenceCopyNumber ?? call.RefPloidy.Value;
                    int CN = call.CN;
                    if (call.AltAllele == "." && optionsSkipDiploid) continue;

                    if (CN > MaxCn) CN = MaxCn;
                    string chr = call.Chr;
                    int overlapStart = Math.Max(call.Start, interval.Start);
                    int overlapEnd = Math.Min(call.End, interval.End);
                    if (overlapStart >= overlapEnd) continue;
                    int overlapBases = overlapEnd - overlapStart;
                    int thisCallBasesExcluded = 0;
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
                            thisCallBasesExcluded += excludeOverlapEnd - excludeOverlapStart;
                            overlapBases -= excludeOverlapEnd - excludeOverlapStart;
                            // if majority of the region is in exclude intervals, don't consider any overlap
                            // N.B.: the denominator here looks dubious -- why compare overlap bases to the excluded interval overlap size,
                            // rather than, say, the size of the original truth interval?  Or the length of the overlap of the current CNV?
                            if (overlapBases / Math.Max(excludeOverlapEnd - excludeOverlapStart, 1) < 0.1)
                            {
                                thisCallBasesExcluded += overlapBases;
                                excludeIntervalBases += overlapBases;
                                overlapBases = 0;
                                break;
                            }
                        }
                    }

                    totalIntervalRefPloidy.Add((refPloidy, overlapBases));

                    if (call.PassFilter || !includePassingOnly)
                    {
                        totalOverlapBases += overlapBases;
                        baseCounter.BaseCount[knownCn, CN, refPloidy] += overlapBases;

                        if (knownCn == CN)
                        {
                            if (CN == refPloidy)
                                thisIntervalBasesTrueNegative += overlapBases;
                            else
                                thisIntervalBasesTruePositive += overlapBases;
                        }
                        else
                        {
                            if (knownCn == refPloidy)
                                thisIntervalBasesFalsePositive += overlapBases;
                            else
                                thisIntervalBasesFalseNegative += overlapBases;
                        }
                        thisIntervalBasesNoCall -= overlapBases;
                        thisIntervalBasesNoCall -= thisCallBasesExcluded;
                        thisIntervalBasesExcluded += thisCallBasesExcluded;
                    }

                    interval.BasesCovered += overlapBases;

                    if (knownCn == CN)
                        interval.BasesCalledCorrectly += overlapBases;
                    else
                        interval.BasesCalledIncorrectly += overlapBases;

                    if (_cnvChecker.RegionsOfInterest.Empty() ||
                        !_cnvChecker.RegionsOfInterest.ContainsKey(chr)) continue;

                    foreach (CNInterval roiInterval in _cnvChecker.RegionsOfInterest[chr])
                    {
                        int roiOverlapStart = Math.Max(roiInterval.Start, overlapStart);
                        int roiOverlapEnd = Math.Min(roiInterval.End, overlapEnd);
                        if (roiOverlapStart >= roiOverlapEnd) continue;
                        int roiOverlapBases = roiOverlapEnd - roiOverlapStart;
                        if (call.PassFilter || !includePassingOnly)
                        {
                            totalRoiOverlapBases += roiOverlapBases;
                            baseCounter.RoiBaseCount[knownCn, CN, refPloidy] += roiOverlapBases;
                        }
                    }
                }

                if (baseCounter.MinSize == 0 && baseCounter.MaxSize > 100000)
                    Console.WriteLine($"Truth {chromosome}:{interval.Start}-{interval.End} CN={knownCn} base counts TP/TN/FP/FN/NC/EXCL {thisIntervalBasesTruePositive} {thisIntervalBasesTrueNegative} {thisIntervalBasesFalsePositive} {thisIntervalBasesFalseNegative} {thisIntervalBasesNoCall} {thisIntervalBasesExcluded}");

                nonOverlapBases -= (totalOverlapBases + excludeIntervalBases);

                if (!interval.ReferenceCopyNumber.HasValue)
                {
                    if (totalIntervalRefPloidy.Empty())
                    {
                        throw new ArgumentException(
                            $"Error: Truth variant {interval.Chromosome}:{interval.Start}-{interval.End} with no overlapping " +
                            $"Canvas calls. Reference ploidy cannot be determined! Please provide reference ploidy via command line options");
                    }

                    interval.ReferenceCopyNumber = Convert.ToInt32(Math.Round(Utilities.WeightedMean(
                        totalIntervalRefPloidy.Select(x => (double)x.ploidy).ToList(),
                        totalIntervalRefPloidy.Select(x => (double)Math.Max(x.length, 1)).ToList())));
                }
                if (nonOverlapBases < 0)
                {
                    throw new InvalidDataException($"Truth variant {interval.Chromosome}:{interval.Start}-{interval.End} has negative non-overlap bases");
                }
                baseCounter.NoCalls[knownCn, interval.ReferenceCopyNumber.Value] += nonOverlapBases;
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