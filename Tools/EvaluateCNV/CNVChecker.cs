using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using CanvasCommon;
using Illumina.Common;
using Isas.SequencingFiles;
using Isas.SequencingFiles.Vcf;

namespace EvaluateCNV
{

    class CNInterval
    {
        public string Chromosome { get; }
        public int Start; // 0-based inclusive
        public int End; // 0-based exclusive
        public int CN;
        public int ReferenceCopyNumber = 2; // updated based on ploidy bed
        public int BasesCovered;
        public int BasesExcluded;
        public int BasesCalledCorrectly;

        public int BasesNotCalled => Length - BasesExcluded - BasesCalledCorrectly - BasesCalledIncorrectly;

        public int Length
        {
            get { return End - Start; }
        }

        public int BasesCalledIncorrectly;

        public override string ToString()
        {
            return $"{Chromosome}:{Start + 1}-{End}";
        }

        public CNInterval(string chromosome)
        {
            Chromosome = chromosome;
        }
    }

    class CNVCall
    {
        public string Chr;
        public int Start; // 0-based inclusive
        public int End; // 0-based exclusive
        public int CN;

        public int Length
        {
            get { return End - Start; }
        }

        public CNVCall(string chr, int start, int end, int cn)
        {
            Chr = chr;
            Start = start;
            End = end;
            CN = cn;
        }

        public int Overlap(PloidyInterval ploidyInterval)
        {
            int overlapStart = Math.Max(Start, ploidyInterval.Start);
            int overlapEnd = Math.Min(End, ploidyInterval.End);
            if (overlapStart >= overlapEnd) return 0;
            return overlapEnd - overlapStart;
        }
    }

    class CNVChecker
    {
        #region Members
        Dictionary<string, List<CNInterval>> KnownCN = null;
        Dictionary<string, List<CNInterval>> RegionsOfInterest = null;
        Dictionary<string, List<CNInterval>> ExcludeIntervals = null;
        public double? DQscoreThreshold { get; }

        public CNVChecker(double? dQscoreThreshold)
        {
            DQscoreThreshold = dQscoreThreshold;
        }
        #endregion

        /// <summary>
        /// Load known CN data from a .bed file.  File lines have fields:
        /// chromosome, start, end, chromcountA, chromcountB
        /// So, copy number is the sum of the last 2 fields, major chromosome count is the max of the last 2 fields.
        /// </summary>
        /// <param name="oracleBedPath"></param>
        protected Dictionary<string, List<CNInterval>> LoadIntervalsFromBed(string oracleBedPath, bool getCN, double heterogeneityFraction)
        {
            bool stripChr = false;
            int count = 0;
            long totalBases = 0;
            Dictionary<string, List<CNInterval>> bedIntervals = new Dictionary<string, List<CNInterval>>();
            using (StreamReader reader = new StreamReader(oracleBedPath))
            {
                while (true)
                {
                    string fileLine = reader.ReadLine();
                    if (fileLine == null) break;
                    if (fileLine.Length == 0 || fileLine[0] == '#') continue;
                    string[] bits = fileLine.TrimEnd('\t').Split('\t');
                    if (bits.Length < 3) continue;
                    string chromosome = bits[0];
                    if (stripChr) chromosome = chromosome.Replace("chr", "");
                    if (!bedIntervals.ContainsKey(chromosome)) bedIntervals[chromosome] = new List<CNInterval>();
                    CNInterval interval = new CNInterval(chromosome);
                    interval.Start = int.Parse(bits[1]);
                    interval.End = int.Parse(bits[2]);
                    if (getCN) // bits.Length >= 5)
                    {
                        if (heterogeneityFraction < 1 && bits.Length > 5 && int.Parse(bits[3]) == 1 && int.Parse(bits[4]) == 1)
                            if (heterogeneityFraction > double.Parse(bits[5]))
                                continue;
                        interval.CN = int.Parse(bits[3]) + int.Parse(bits[4]);
                    }
                    totalBases += interval.Length;
                    bedIntervals[chromosome].Add(interval);
                    count++;
                }
            }
            Console.WriteLine(">>>Loaded {0} CN intervals ({1} bases)", count, totalBases);
            return bedIntervals;
        }


        protected void LoadKnownCNVCF(string oracleVCFPath)
        {
            bool stripChr = false;

            // Load our "oracle" of known copy numbers:
            this.KnownCN = new Dictionary<string, List<CNInterval>>();
            int count = 0;
            using (GzipReader reader = new GzipReader(oracleVCFPath))
            {
                while (true)
                {
                    string fileLine = reader.ReadLine();
                    if (fileLine == null) break;
                    if (fileLine.Length == 0 || fileLine[0] == '#') continue;
                    string[] bits = fileLine.Split('\t');
                    string chromosome = bits[0];
                    if (stripChr) chromosome = chromosome.Replace("chr", "");
                    if (!KnownCN.ContainsKey(chromosome)) KnownCN[chromosome] = new List<CNInterval>();
                    CNInterval interval = new CNInterval(chromosome);
                    interval.Start = int.Parse(bits[1]);
                    interval.CN = -1;
                    string[] infoBits = bits[7].Split(';');
                    foreach (string subBit in infoBits)
                    {
                        if (subBit.StartsWith("CN="))
                        {
                            float tempCN = float.Parse(subBit.Substring(3));
                            if (subBit.EndsWith(".5"))
                            {
                                interval.CN = (int)Math.Round(tempCN + 0.1); // round X.5 up to X+1
                            }
                            else
                            {
                                interval.CN = (int)Math.Round(tempCN); // Round off
                            }
                        }
                        if (subBit.StartsWith("END="))
                        {
                            interval.End = int.Parse(subBit.Substring(4));
                        }
                    }
                    // Parse CN from Canvas output:
                    if (bits.Length > 8)
                    {
                        string[] subBits = bits[8].Split(':');
                        string[] subBits2 = bits[9].Split(':');
                        for (int subBitIndex = 0; subBitIndex < subBits.Length; subBitIndex++)
                        {
                            if (subBits[subBitIndex] == "CN")
                            {
                                interval.CN = int.Parse(subBits2[subBitIndex]);
                            }
                        }
                    }
                    if (interval.End == 0 || interval.CN < 0)
                    {
                        Console.WriteLine("Error - bogus record!");
                        Console.WriteLine(fileLine);
                    }
                    else
                    {
                        KnownCN[chromosome].Add(interval);
                        count++;
                    }
                }
            }
            Console.WriteLine(">>>Loaded {0} known-CN intervals", count);
        }

        protected void LoadRegionsOfInterest(string bedPath)
        {
            if (string.IsNullOrEmpty(bedPath)) return;
            if (!File.Exists(bedPath))
            {
                throw new ArgumentException(string.Format("* Error: ROI bed file not found at '{0}'", bedPath));
            }
            this.RegionsOfInterest = this.LoadIntervalsFromBed(bedPath, false, 1.0);
            List<string> keys = this.RegionsOfInterest.Keys.ToList();
            foreach (string key in keys)
            {
                this.RegionsOfInterest[string.Format("chr{0}", key)] = this.RegionsOfInterest[key];
            }
        }

        protected void LoadKnownCN(string oraclePath, double heterogeneityFraction)
        {
            if (!File.Exists(oraclePath))
            {
                throw new ArgumentException(string.Format("* Error: Truth vcf not found at '{0}'", oraclePath));
            }

            if (oraclePath.EndsWith(".bed"))
            {
                this.KnownCN = this.LoadIntervalsFromBed(oraclePath, true, heterogeneityFraction);
                return;
            }
            LoadKnownCNVCF(oraclePath);
            SummarizeTruthSetStatistics();
        }

        protected void SummarizeTruthSetStatistics()
        {
            List<long> EventSizes = new List<long>();
            double MeanEventSize = 0;
            int countUnder10KB = 0;
            int count10kb50kb = 0;
            int count50kb500kb = 0;
            int count500kbplus = 0;
            foreach (string key in KnownCN.Keys)
            {
                foreach (CNInterval interval in KnownCN[key])
                {
                    if (interval.CN == 2) continue;
                    long length = interval.Length;
                    EventSizes.Add(length);
                    MeanEventSize += length;
                    if (length <= 10000)
                    {
                        countUnder10KB++;
                    }
                    else if (length <= 50000)
                    {
                        count10kb50kb++;
                    }
                    else if (length <= 500000)
                    {
                        count50kb500kb++;
                    }
                    else
                    {
                        count500kbplus++;
                    }
                }
            }
            EventSizes.Sort();
            MeanEventSize /= EventSizes.Count;
            Console.WriteLine("Known CN: Mean length {0:F4}", MeanEventSize);
            Console.WriteLine("up to 10kb: {0}", countUnder10KB);
            Console.WriteLine("10kb-50kb: {0}", count10kb50kb);
            Console.WriteLine("50kb-500kb: {0}", count50kb500kb);
            Console.WriteLine("500kb+: {0}", count500kbplus);
            if (EventSizes.Count > 0)
                Console.WriteLine("Median size: {0}", EventSizes[EventSizes.Count / 2]);
        }


        protected int GetCopyNumber(VcfVariant variant, out int end)
        {
            int CN = -1;
            end = -1;
            if (variant.GenotypeColumns != null && variant.GenotypeColumns.Count > 0)
            {
                Dictionary<string, string> genotype = variant.GenotypeColumns[variant.GenotypeColumns.Count - 1];
                if (genotype.ContainsKey("CN"))
                {
                    CN = int.Parse(genotype["CN"]);
                }
                if (genotype.ContainsKey("END"))
                {
                    end = int.Parse(genotype["END"]);
                }
            }
            if (variant.InfoFields.ContainsKey("END"))
            {
                end = int.Parse(variant.InfoFields["END"]);
            }
            if (variant.InfoFields.ContainsKey("CN"))
            {
                CN = int.Parse(variant.InfoFields["CN"]);
            }

            return CN;
        }

        protected void CountExcludedBasesInTruthSetIntervals()
        {
            foreach (string key in KnownCN.Keys)
            {
                if (!ExcludeIntervals.ContainsKey(key)) continue;
                foreach (CNInterval interval in KnownCN[key])
                {
                    foreach (CNInterval excludeInterval in ExcludeIntervals[key])
                    {
                        int overlapStart = Math.Max(interval.Start, excludeInterval.Start);
                        int overlapEnd = Math.Min(interval.End, excludeInterval.End);
                        if (overlapStart >= overlapEnd) continue;
                        interval.BasesExcluded += overlapEnd - overlapStart;
                        //Console.WriteLine("Interval {0}:{1}-{2} excludes {3} bases due to overlap with excluded interval {4}:{5}-{6}",
                        //    key, interval.Start, interval.End, overlapEnd - overlapStart,
                        //    key, excludeInterval.Start, excludeInterval.End);
                    }
                }
            }
        }

        protected IEnumerable<CNVCall> GetCnvCallsFromVcf(string vcfPath, bool includePassingOnly)
        {
            using (VcfReader reader = new VcfReader(vcfPath, false))
            {
                if (DQscoreThreshold.HasValue)
                {
                    var match = reader.HeaderLines.FirstOrDefault(stringToCheck => stringToCheck.Contains("DQSCORE"));
                    if (match == null)
                        throw new ArgumentException($"File {vcfPath} does not contain DQscore INFO field.");
                }

                foreach (VcfVariant variant in reader.GetVariants())
                {

                    int end;
                    int CN = GetCopyNumber(variant, out end);
                    if (includePassingOnly && variant.Filters != "PASS") continue;
                    if (DQscoreThreshold.HasValue && !variant.InfoFields.ContainsKey("DQSCORE") && CN == 2)
                    {
                        yield return new CNVCall(variant.ReferenceName, variant.ReferencePosition, end, CN);
                        continue;
                    }
                    if (DQscoreThreshold.HasValue && variant.InfoFields.ContainsKey("DQSCORE") &&
                        double.Parse(variant.InfoFields["DQSCORE"]) < DQscoreThreshold.Value)
                        continue;
                    yield return new CNVCall(variant.ReferenceName, variant.ReferencePosition, end, CN);
                }
            }
        }

        protected IEnumerable<CNVCall> GetCnvCallsFromBed(string bedPath, int[] cnIndices = null)
        {
            if (cnIndices == null) { cnIndices = new int[] { 3 }; }
            int maxCnIndex = cnIndices.Max();
            using (StreamReader reader = new StreamReader(bedPath))
            {
                string line;
                string[] toks;
                while ((line = reader.ReadLine()) != null)
                {
                    if (line.StartsWith("#")) { continue; } // skip comments
                    toks = line.Split('\t');
                    if (toks.Length <= maxCnIndex)
                    {
                        Console.WriteLine("Error: Line has fewer than {0} columns: {1}", maxCnIndex + 1, line);
                        continue;
                    }
                    string chr;
                    int start;
                    int end;
                    int cn;
                    try
                    {
                        chr = toks[0];
                        start = int.Parse(toks[1]);
                        end = int.Parse(toks[2]);
                        cn = cnIndices.Sum(cnIndex => int.Parse(toks[cnIndex]));
                    }
                    catch
                    {
                        Console.WriteLine("Error: Failed to parse line: {0}", line);
                        continue;
                    }
                    yield return new CNVCall(chr, start, end, cn);
                }
            }
        }

        protected void ComputeAccuracy(string truthSetPath, string cnvCallsPath, StreamWriter outputWriter, PloidyInfo ploidyInfo, bool includePassingOnly)
        {
            int totalVariants = 0;
            long totalVariantBases = 0;
            int maxCN = 10;
            //BaseCount[TrueCN, CalledCN] = total # of bases of overlap for this tuple
            long[,] BaseCount = new long[maxCN + 1, maxCN + 1];
            long[,] ROIBaseCount = new long[maxCN + 1, maxCN + 1];

            // Make a note of how many bases in the truth set are not *actually* considered to be known bases, using
            // the "cnaqc" exclusion set:
            this.CountExcludedBasesInTruthSetIntervals();
            if (DQscoreThreshold.HasValue && !Path.GetFileName(cnvCallsPath).ToLower().Contains("vcf"))
                throw new ArgumentException("CNV.vcf must be in a vcf format when --dqscore option is used");
            IEnumerable<CNVCall> calls = Path.GetFileName(cnvCallsPath).ToLower().Contains("vcf")
                ? GetCnvCallsFromVcf(cnvCallsPath, includePassingOnly)
                : GetCnvCallsFromBed(cnvCallsPath);

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
                if (variantBasesOverlappingPloidyRegion > 0 || (CN != 2 && variantBasesOverlappingPloidyRegion < call.Length))
                    totalVariants++;

                if (CN > maxCN) CN = maxCN;
                string chr = call.Chr;
                if (!KnownCN.ContainsKey(chr)) chr = call.Chr.Replace("chr", "");
                if (!KnownCN.ContainsKey(chr)) chr = "chr" + call.Chr;
                if (!KnownCN.ContainsKey(chr))
                {
                    Console.WriteLine("Error: Skipping variant call for chromosome {0} with no truth data", call.Chr);
                    continue;
                }
                foreach (CNInterval interval in KnownCN[chr])
                {
                    int overlapStart = Math.Max(call.Start, interval.Start);
                    int overlapEnd = Math.Min(call.End, interval.End);
                    if (overlapStart >= overlapEnd) continue;
                    int overlapBases = overlapEnd - overlapStart;
                    // We've got an overlap interval.  Kill off some bases from this interval, if it happens
                    // to overlap with an excluded interval:
                    if (ExcludeIntervals.ContainsKey(chr))
                    {
                        foreach (CNInterval excludeInterval in ExcludeIntervals[chr])
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
                    interval.BasesCovered += overlapBases;
                    if (knownCN == CN)
                    {
                        interval.BasesCalledCorrectly += overlapBases;
                    }
                    else
                    {
                        interval.BasesCalledIncorrectly += overlapBases;
                    }

                    if (this.RegionsOfInterest != null && this.RegionsOfInterest.ContainsKey(chr))
                    {
                        foreach (CNInterval roiInterval in this.RegionsOfInterest[chr])
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
            foreach (string chr in KnownCN.Keys)
            {
                foreach (CNInterval interval in KnownCN[chr])
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

            IEnumerable<CNInterval> allIntervals = KnownCN.SelectMany(kvp => kvp.Value);

            // find truth interval with highest number of false negatives (hurts recall)
            IEnumerable<CNInterval> variantIntervals = allIntervals.Where(interval => interval.CN != interval.ReferenceCopyNumber);
            CNInterval intervalMaxFalseNegatives = variantIntervals.MaxBy(interval => interval.BasesNotCalled + interval.BasesCalledIncorrectly);
            Console.WriteLine($"Truth interval with most false negatives (hurts recall): {intervalMaxFalseNegatives}");

            // find truth interval with highest number of false positive (hurts precision)
            IEnumerable<CNInterval> refIntervals = allIntervals.Where(interval => interval.CN == interval.ReferenceCopyNumber);
            CNInterval intervalMaxFalsePositives = refIntervals.MaxBy(interval => interval.BasesCalledIncorrectly);
            Console.WriteLine($"Truth interval with most false positives (hurts precision): {intervalMaxFalsePositives}");

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
                    long bases = BaseCount[trueCN, callCN];
                    totalBases += bases;
                    if (trueCN == callCN) totalBasesRight += bases;
                    if ((trueCN < 2 && callCN < 2) || (trueCN == 2 && callCN == 2) ||
                        (trueCN > 2 && callCN > 2))
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
            for (int trueCN = 0; trueCN <= maxCN; trueCN++)
            {
                for (int callCN = 0; callCN <= maxCN; callCN++)
                {
                    long bases = ROIBaseCount[trueCN, callCN];
                    ROIBases += bases;
                    if (trueCN == callCN) ROIBasesCorrect += bases;
                    if ((trueCN < 2 && callCN < 2) || (trueCN == 2 && callCN == 2) || (trueCN > 2 && callCN > 2)) ROIBasesCorrectDirection += bases;
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
            outputWriter.WriteLine("Accuracy\t{0:F4}", 100 * totalBasesRight / (double)totalBases);
            outputWriter.WriteLine("DirectionAccuracy\t{0:F4}", 100 * totalBasesRightDirection / (double)totalBases);
            outputWriter.WriteLine("Recall\t{0:F4}",
                100 * (isGainBasesCorrect + isLossBasesCorrect) / (double)(isGainBases + isLossBases));
            outputWriter.WriteLine("DirectionRecall\t{0:F4}",
                100 * (isGainBasesCorrectDirection + isLossBasesCorrectDirection) / (double)(isGainBases + isLossBases));
            outputWriter.WriteLine("Precision\t{0:F4}",
                100 * (isGainBasesCorrect + isLossBasesCorrect) / (double)(callGainBases + callLossBases));
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
            outputWriter.WriteLine("MeanEventAccuracy\t{0:F4}", 100 * meanAccuracy);
            outputWriter.WriteLine("MedianEventAccuracy\t{0:F4}", 100 * medianAccuracy);
            outputWriter.WriteLine("VariantEventsCalled\t{0}", totalVariants);
            outputWriter.WriteLine("VariantBasesCalled\t{0}", totalVariantBases);
            if (ROIBases > 0)
            {
                outputWriter.WriteLine("ROIAccuracy\t{0:F4}", 100 * ROIBasesCorrect / (double)ROIBases);
                outputWriter.WriteLine("ROIDirectionAccuracy\t{0:F4}", 100 * ROIBasesCorrectDirection / (double)ROIBases);
            }
            // to separate passing and all variant results
            outputWriter.WriteLine();
        }

        public void Evaluate(string truthSetPath, string cnvCallsPath, string excludedBed, string outputPath, EvaluateCnvOptions options)
        {
            double heterogeneityFraction = options.HeterogeneityFraction;
            var ploidyInfo = PloidyInfo.LoadPloidyFromBedFile(options.PloidyBed?.FullName);

            LoadKnownCN(truthSetPath, heterogeneityFraction);
            ploidyInfo.MakeChromsomeNameAgnosticWithAllChromosomes(KnownCN.Keys);
            SetTruthsetReferencePloidy(ploidyInfo);

            // LoadRegionsOfInterest(options.RoiBed?.FullName);
            if (!string.IsNullOrEmpty(excludedBed))
            {
                this.ExcludeIntervals = LoadIntervalsFromBed(excludedBed, false, 1.0);
                // cheesy logic to handle different chromosome names:
                List<string> keys = this.ExcludeIntervals.Keys.ToList();
                foreach (string key in keys)
                {
                    this.ExcludeIntervals[key.Replace("chr", "")] = this.ExcludeIntervals[key];
                }
            }
            Console.WriteLine("TruthSet\t{0}", truthSetPath);
            Console.WriteLine("CNVCalls\t{0}", cnvCallsPath);
            using (StreamWriter outputWriter = new StreamWriter(outputPath))
            {
                if (Path.GetFileName(cnvCallsPath).ToLower().Contains("vcf"))
                {
                    ComputeAccuracy(truthSetPath, cnvCallsPath, outputWriter, ploidyInfo, true);
                }
                ComputeAccuracy(truthSetPath, cnvCallsPath, outputWriter, ploidyInfo, false);
            }
            Console.WriteLine(">>>Done - results written to {0}", outputPath);
        }

        private void SetTruthsetReferencePloidy(PloidyInfo ploidyInfo)
        {
            foreach (string chromosome in KnownCN.Keys)
            {
                foreach (CNInterval truthInterval in KnownCN[chromosome])
                {
                    foreach (PloidyInterval ploidyRegion in ploidyInfo.PloidyByChromosome[chromosome])
                    {
                        // truth interval must be completely contained within the ploidy region
                        if (truthInterval.End >= ploidyRegion.Start && truthInterval.Start <= ploidyRegion.End)
                        {
                            truthInterval.ReferenceCopyNumber = ploidyRegion.Ploidy;
                            break;
                        }
                        if ((truthInterval.Start >= ploidyRegion.Start && truthInterval.Start <= ploidyRegion.End) ||
                            (truthInterval.End >= ploidyRegion.Start && truthInterval.End <= ploidyRegion.End))
                            throw new ApplicationException($"Truth interval {truthInterval} crosses reference ploidy region {ploidyRegion}. Update truth interval");
                    }
                }
            }
        }
    }
}
