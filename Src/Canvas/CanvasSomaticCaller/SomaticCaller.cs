using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using SequencingFiles;
using CanvasCommon;
using SequencingFiles.Vcf;

namespace CanvasSomaticCaller
{
    public class SomaticCaller
    {

        #region Members
        // Static:

        // Data:
        List<CanvasSegment> Segments;
        List<SegmentPloidy> AllPloidies;
        CopyNumberOracle CNOracle = null;
        CoveragePurityModel Model;
        Dictionary<string, List<GenomicBin>> ExcludedIntervals = new Dictionary<string, List<GenomicBin>>();
        List<long> HeterogeneousSegmentsSignature = new List<long>();

        // File paths:
        public string TruthDataPath;
        public string SomaticVCFPath;
        protected string OutputFolder;
        protected string TempFolder;

        // Weighting parameters:
        // derived from fitting logistic regression to "true" and "false" ploidy and purity models
        // on tumour sequenced dataset of ~50 samples
        //public double PercentNormal2WeightingFactor = 0.3;
        //public double DeviationScoreWeightingFactor = 0.25;
        //public double CN2WeightingFactor = 0.35;
        //public double DiploidDistanceScoreWeightingFactor = 0.25;
        // New weights - incrementally better separation between top-scoring true model
        // and the next runner-up:

        // Parameters:
        public float? userPloidy;
        public float? userPurity;
        protected float MeanCoverage = 30;
        private double CoverageWeightingFactor; // Computed from CoverageWeighting
        public bool IsEnrichment = false;
        public bool IsDbsnpVcf = false;
        public bool IsTrainingMode = false;
        protected PloidyInfo ReferencePloidy;
        public SomaticCallerParameters somaticCallerParameters;

        public bool FFPEMode = false; // Assume MAF and Coverage are independent/uncorrelated in FFPEMode (always false for now)
        private const double EMPosteriorProbThres = 0.01; // Controls whether a segment contributes to Mu and Sigma estimates
        private const double EMOmegaThres = 0.01; // Controls when to update means
        private const double EMLikelihoodThres = 1; // Controls when to update means
        public int QualityFilterThreshold { get; set; } = 10;

        #endregion

        /// <summary>
        /// Load the expected ploidy for sex chromosomes from a .bed file.  This lets us know that, for instance, copy number 2
        /// on chrX is a GAIN (not REF) call for a male (XY) sample.
        /// </summary>
        public void LoadReferencePloidy(string filePath)
        {
            Console.WriteLine(">>>LoadReferencePloidy({0})", filePath);
            ReferencePloidy = PloidyInfo.LoadPloidyFromBedFile(filePath);
        }

        /// <summary>
        /// Setup: Model various copy ploidies.
        /// </summary>
        private void InitializePloidies()
        {
            Console.WriteLine("{0} Initialize ploidy models...", DateTime.Now);
            this.AllPloidies = new List<SegmentPloidy>();
            double diploidPredictedMAF = CanvasCommon.Utilities.EstimateDiploidMAF(2, this.MeanCoverage);
            for (int copyNumber = 0; copyNumber <= somaticCallerParameters.MaximumCopyNumber; copyNumber++)
            {
                for (int majorCount = copyNumber; majorCount * 2 >= copyNumber; majorCount--)
                {
                    SegmentPloidy ploidy = new SegmentPloidy();
                    ploidy.CopyNumber = copyNumber;
                    ploidy.MajorChromosomeCount = majorCount;
                    ploidy.ID = AllPloidies.Count;
                    AllPloidies.Add(ploidy);

                    if (copyNumber == 0)
                    {
                        ploidy.MinorAlleleFrequency = CanvasCommon.Utilities.EstimateDiploidMAF(1, this.MeanCoverage);
                        continue;
                    }
                    float VF = majorCount / (float)copyNumber;
                    ploidy.MinorAlleleFrequency = (VF < 0.5 ? VF : 1 - VF);
                    if (majorCount * 2 == copyNumber)
                    {
                        ploidy.MinorAlleleFrequency = CanvasCommon.Utilities.EstimateDiploidMAF(copyNumber, this.MeanCoverage);
                    }

                }
            }
            Console.WriteLine("{0} Ploidy models prepared.", DateTime.Now);
        }

        /// <summary>
        /// Developer debug option:
        /// Given known CN data, determine the ploidy for each segment, and report this in such a way that we 
        /// can review scatterplots of (MAF, Coverage) for each ploidy.  
        /// </summary>
        private void DebugModelSegmentsByPloidy()
        {
            Console.WriteLine("{0} DebugModelSegmentsByPloidy...", DateTime.Now);
            using (StreamWriter writer = new StreamWriter(Path.Combine(this.OutputFolder, "SegmentsByPloidy.txt")))
            {
                writer.WriteLine("#Chr\tStart\tEnd\tCN\tPloidy\tMAF\tCoverage");
                foreach (CanvasSegment segment in Segments)
                {
                    if (segment.End - segment.Begin < 5000) continue; // skip over itty bitty segments.
                    int CN = this.GetKnownCNForSegment(segment);
                    if (CN < 0) continue;
                    List<float> MAF = new List<float>();
                    foreach (float VF in segment.VariantFrequencies)
                    {
                        MAF.Add(VF > 0.5 ? 1 - VF : VF);
                    }
                    if (MAF.Count < somaticCallerParameters.MinimumVariantFrequenciesForInformativeSegment) continue;
                    MAF.Sort();
                    float MedianMAF = MAF[MAF.Count / 2];
                    int medianCoverage = (int)Math.Round(CanvasCommon.Utilities.Median(segment.Counts));
                    // Identify the most plausible ploidy:
                    string ploidy = "?";
                    switch (CN)
                    {
                        case 0:
                            ploidy = "-";
                            break;
                        case 1:
                            ploidy = "A";
                            break;
                        case 2:
                            ploidy = "AB";
                            if (MedianMAF < 0.25) ploidy = "AA";
                            break;
                        case 3:
                            ploidy = "ABB";
                            if (MedianMAF < 0.165) ploidy = "AAA";
                            break;
                        case 4:
                            ploidy = "AABB";
                            if (MedianMAF < 0.375) ploidy = "AAAB";
                            if (MedianMAF < 0.125) ploidy = "AAAA";
                            break;
                    }
                    writer.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t", segment.Chr, segment.Begin, segment.End,
                        CN, ploidy, MedianMAF, medianCoverage);
                }
            }
            Console.WriteLine("{0} DebugModelSegmentsByPloidy complete", DateTime.Now);
        }

        private void DebugModelSegmentCoverageByCN()
        {
            Console.WriteLine("{0} DebugModelSegmentCoverageByCN...", DateTime.Now);
            int histogramBinCount = 1024;
            // Initialize histograms:
            int[][] CNHistogram = new int[10][];
            for (int CN = 0; CN < 10; CN++)
            {
                CNHistogram[CN] = new int[histogramBinCount];
            }
            int[] histogramAll = new int[histogramBinCount];

            // Accumulate histograms:
            foreach (CanvasSegment segment in Segments)
            {
                if (segment.End - segment.Begin < 5000) continue; // skip over itty bitty segments.
                int count = (int)Math.Round(CanvasCommon.Utilities.Median(segment.Counts));
                if (count >= histogramBinCount)
                {
                    continue;
                }
                else if (count < 0)
                {
                    continue;
                }
                int CN = this.GetKnownCNForSegment(segment);
                histogramAll[count]++;
                if (CN < 0 || CN >= CNHistogram.Length) continue;
                CNHistogram[CN][count]++;
            }

            // Report summary stats:
            for (int CN = 0; CN < 10; CN++)
            {
                float total = 0;
                foreach (int value in CNHistogram[CN]) total += value;
                if (total == 0) continue;
                int runningTotal = 0;
                int median = -1;
                float mean = 0;
                for (int bin = 0; bin < CNHistogram[CN].Length; bin++)
                {
                    runningTotal += CNHistogram[CN][bin];
                    if (median < 0 && runningTotal * 2 >= total)
                    {
                        median = bin;
                    }
                    mean += bin * CNHistogram[CN][bin];
                }
                mean /= total;
                double stddev = 0;
                for (int bin = 0; bin < CNHistogram[CN].Length; bin++)
                {
                    float diff = bin - mean;
                    stddev += (diff * diff) * CNHistogram[CN][bin];
                }
                stddev = Math.Sqrt(stddev / total);
                Console.WriteLine("CN {0} median {1} mean {2:F2} stddev {3:F2}", CN, median, mean, stddev);
            }
            // Dump histograms to a text file:
            using (StreamWriter writer = new StreamWriter(Path.Combine(this.OutputFolder, "SegmentCoverageByCN.txt")))
            {
                writer.Write("#Bin\tAllSegments\t");
                for (int CN = 0; CN < 10; CN++)
                {
                    writer.Write("CN{0}\t", CN);
                }
                writer.WriteLine();
                for (int bin = 0; bin < histogramBinCount; bin++)
                {
                    writer.Write("{0}\t", bin);
                    writer.Write("{0}\t", histogramAll[bin]);
                    for (int CN = 0; CN < 10; CN++)
                    {
                        writer.Write("{0}\t", CNHistogram[CN][bin]);
                    }
                    writer.WriteLine();
                }
            }
            Console.WriteLine("{0} DebugModelSegmentCoverageByCN complete", DateTime.Now);
        }

        /// <summary>
        /// Developer troubleshooting method:
        /// Given truth data about CNs, let's derive and log a histogram of bin-coverage for each copy number.
        /// We'll also review summary statistics by CN.
        /// </summary>
        /// <returns></returns>
        private void DebugModelCoverageByCN()
        {
            Console.WriteLine("{0} DebugModelCoverageByCN...", DateTime.Now);
            int histogramBinCount = 1024;

            // Initialize histograms:
            int[][] CNHistogram = new int[10][];
            for (int CN = 0; CN < 10; CN++)
            {
                CNHistogram[CN] = new int[histogramBinCount];
            }
            int[] histogramAll = new int[histogramBinCount];

            // Accumulate histograms:
            foreach (CanvasSegment segment in Segments)
            {
                foreach (float tempCount in segment.Counts)
                {
                    int count = (int)Math.Round(tempCount);
                    if (count >= histogramBinCount)
                    {
                        continue;
                    }
                    else if (count < 0)
                    {
                        continue;
                    }
                    histogramAll[count]++;
                }
                int CN = this.GetKnownCNForSegment(segment);
                if (CN < 0 || CN >= CNHistogram.Length) continue;
                foreach (float tempCount in segment.Counts)
                {
                    int count = (int)Math.Round(tempCount);
                    if (count >= histogramBinCount)
                    {
                        //Console.Error.WriteLine("Skipping count {0} >= histogram bin count {1}.", count, histogramBinCount);
                        continue;
                    }
                    else if (count < 0)
                    {
                        //Console.Error.WriteLine("Skipping count {0} < 0.", count);
                        continue;
                    }
                    CNHistogram[CN][count]++;
                }
            }

            // Report summary stats:
            for (int CN = 0; CN < 10; CN++)
            {
                float total = 0;
                foreach (int value in CNHistogram[CN]) total += value;
                if (total == 0) continue;
                int runningTotal = 0;
                int median = -1;
                float mean = 0;
                for (int bin = 0; bin < CNHistogram[CN].Length; bin++)
                {
                    runningTotal += CNHistogram[CN][bin];
                    if (median < 0 && runningTotal * 2 >= total)
                    {
                        median = bin;
                    }
                    mean += bin * CNHistogram[CN][bin];
                }
                mean /= total;
                double stddev = 0;
                for (int bin = 0; bin < CNHistogram[CN].Length; bin++)
                {
                    float diff = bin - mean;
                    stddev += (diff * diff) * CNHistogram[CN][bin];
                }
                stddev = Math.Sqrt(stddev / total);
                Console.WriteLine("CN {0} median {1} mean {2:F2} stddev {3:F2}", CN, median, mean, stddev);
            }
            // Dump histograms to a text file:
            using (StreamWriter writer = new StreamWriter(Path.Combine(this.TempFolder, "BinCoverageByCN.txt")))
            {
                writer.Write("#Bin\tAllBins\t");
                for (int CN = 0; CN < 10; CN++)
                {
                    writer.Write("CN{0}\t", CN);
                }
                writer.WriteLine();
                for (int bin = 0; bin < histogramBinCount; bin++)
                {
                    writer.Write("{0}\t", bin);
                    writer.Write("{0}\t", histogramAll[bin]);
                    for (int CN = 0; CN < 10; CN++)
                    {
                        writer.Write("{0}\t", CNHistogram[CN][bin]);
                    }
                    writer.WriteLine();
                }
            }
            Console.WriteLine("{0} DebugModelCoverageByCN complete", DateTime.Now);
        }

        public void LoadBedFile(string bedPath)
        {
            if (string.IsNullOrEmpty(bedPath)) return;
            this.ExcludedIntervals = CanvasCommon.Utilities.LoadBedFile(bedPath);
        }

        public int CallVariants(string inFile, string variantFrequencyFile, string outputVCFPath, string referenceFolder, string name, double? localSDmertic, CanvasSomaticClusteringMode clusteringMode)
        {
            this.OutputFolder = Path.GetDirectoryName(outputVCFPath);
            this.TempFolder = Path.GetDirectoryName(inFile);
            Console.WriteLine("{0} CallVariants start:", DateTime.Now);
            this.Segments = CanvasSegment.ReadSegments(inFile);

            // Special logic: Increase the allowed model deviation for targeted data.
            if (this.Segments.Count < 500)
                somaticCallerParameters.DeviationFactor = 2.0f;

            // Some debugging output, for developer usage:
            if (!string.IsNullOrEmpty(this.TruthDataPath))
            {
                this.CNOracle = new CopyNumberOracle();
                this.CNOracle.LoadKnownCN(this.TruthDataPath);
            }
            if (this.CNOracle != null)
            {
                this.DebugModelCoverageByCN();
                this.DebugModelSegmentCoverageByCN();
            }

            this.MeanCoverage = CanvasIO.LoadVariantFrequencies(variantFrequencyFile, this.Segments);
            if (this.IsDbsnpVcf)
            {
                int tmpMinimumVariantFreq = somaticCallerParameters.MinimumVariantFrequenciesForInformativeSegment;
                CanvasCommon.Utilities.PruneVariantFrequencies(this.Segments, this.TempFolder, ref tmpMinimumVariantFreq);
                somaticCallerParameters.MinimumVariantFrequenciesForInformativeSegment = tmpMinimumVariantFreq;
            }

            this.InitializePloidies();

            if (this.CNOracle != null) this.DebugModelSegmentsByPloidy();
            List<string> ExtraHeaders = new List<string>();
            try
            {
                ExtraHeaders = CallCNVUsingSNVFrequency(localSDmertic, referenceFolder, clusteringMode);
            }
            catch (Exception e)
            {
                // In a training mode (INTERNAL) somatic model is initialized with a large number of parameter trials. 
                // Some of them might lead to exception as they would fall outside testable range. 
                // For such cases when the IsTrainingMode is set, the program will terminate normally but will produce an empty vcf file. 
                // This will penalize a parameter combination that lead to exception thereby preventing it from creeping into default SomaticCallerParameters.json.
                if (this.IsTrainingMode)
                {
                    Console.WriteLine("Not calling any CNVs. Reason: {0}", e.Message);
                    Segments.Clear();
                    CanvasSegment.WriteSegments(outputVCFPath, this.Segments, Model.DiploidCoverage, referenceFolder, name, ExtraHeaders,
                        this.ReferencePloidy, QualityFilterThreshold);
                    Environment.Exit(0);
                }
                else
                {
                    if (e is SomaticCaller.UncallableDataException)
                    {
                        Console.Error.WriteLine("Not calling any CNVs. Reason: {0}", e.Message);
                        Segments.Clear();
                    }
                    else
                    {
                        throw;
                    }
                }
            }

            string coverageOutputPath = CanvasCommon.Utilities.GetCoverageAndVariantFrequencyOutputPath(outputVCFPath);
            CanvasSegment.WriteCoveragePlotData(this.Segments, this.Model?.DiploidCoverage, this.ReferencePloidy, coverageOutputPath, referenceFolder);

            if (this.ReferencePloidy != null && !string.IsNullOrEmpty(this.ReferencePloidy.HeaderLine))
            {
                ExtraHeaders.Add(this.ReferencePloidy.HeaderLine);
            }

            CanvasSegment.AssignQualityScores(this.Segments, CanvasSegment.QScoreMethod.Logistic);
            this.FilterCNVCalls();

            // Merge *neighboring* segments that got the same copy number call.
            // Enrichment is not allowed to merge non-adjacent segments, since many of those merges would
            // jump across non-manifest intervals.
            if (this.IsEnrichment)
            {
                CanvasSegment.MergeSegments(ref this.Segments, somaticCallerParameters.MinimumCallSize, 1);
            }
            else
            {
                CanvasSegment.MergeSegmentsUsingExcludedIntervals(ref this.Segments, somaticCallerParameters.MinimumCallSize, ExcludedIntervals);
            }

            if (this.CNOracle != null)
            {
                this.DebugEvaluateCopyNumberCallAccuracy();
                this.GenerateReportVersusKnownCN();
                this.GenerateExtendedReportVersusKnownCN();
            }

            ExtraHeaders.Add($"##EstimatedChromosomeCount={this.EstimateChromosomeCount():F2}");

            // Write out results:
            CanvasSegment.WriteSegments(outputVCFPath, this.Segments, Model.DiploidCoverage, referenceFolder, name, ExtraHeaders, this.ReferencePloidy, QualityFilterThreshold);

            return 0;
        }

        /// <summary>
        /// Set segment.Filter for each of our segments.
        /// </summary>
        protected void FilterCNVCalls()
        {
            string qualityFilter = $"q{QualityFilterThreshold}";
            foreach (var segment in this.Segments)
            {
                string filter = null;
                if (segment.QScore < QualityFilterThreshold)
                {
                    filter = qualityFilter;
                }
                if (segment.End - segment.Begin < 10000)
                {
                    if (filter != null)
                        filter = filter + ";L10kb";
                    else
                        filter = "L10kb";
                }
                if (filter == null)
                    filter = "PASS";

                segment.Filter = filter;
            }
        }

        /// <summary>
        /// Developer option:
        /// For each segment with known CN and a reasonably large number of variants, identify the median minor allele
        /// frequency.  Using that information, find the ploidy + purity that best fits the data.  
        /// so skip over those.  For all informative segments, note the purity; track whether we see a consistent 
        /// purity level.
        /// </summary>
        protected void DerivePurityEstimateFromVF()
        {
            double totalPurity = 0;
            double totalWeight = 0;
            int attemptedSegmentCount = 0;
            int informativeSegmentCount = 0;
            List<double> allPurities = new List<double>();
            foreach (CanvasSegment segment in this.Segments)
            {
                int CN = this.GetKnownCNForSegment(segment);
                // Require the segment have a known CN and reasonably large number of variants:
                if (CN < 0) continue;
                if (segment.VariantFrequencies.Count < somaticCallerParameters.MinimumVariantFrequenciesForInformativeSegment) continue;

                List<float> MAF = new List<float>();
                foreach (float VF in segment.VariantFrequencies)
                {
                    MAF.Add(VF > 0.5 ? 1 - VF : VF);
                }
                MAF.Sort();
                float MedianMAF = MAF[MAF.Count / 2];

                // Note that segments with CN=2 or CN=4 and MAF ~= 50% aren't informative here (purity could be anywhere 
                // from 0% and 100% and we wouldn't know the difference).  Skip non-informative segments:
                if (CN == 2 && Math.Abs(0.5 - MedianMAF) < 0.15) continue;
                if (CN == 4 && Math.Abs(0.5 - MedianMAF) < 0.1) continue;
                attemptedSegmentCount++;
                // Skip over segments that have unexpectedly broad distributions of MAF - in some cases we see
                // what look like bimodal histograms.  (Maybe these are segments that don't have consistent CN, or
                // that have consistent CN but inconsistent ploidy)
                int countBad = 0;
                int countGood = 0;
                foreach (float value in MAF)
                {
                    if (Math.Abs(MedianMAF - value) < 0.05)
                    {
                        countGood++;
                    }
                    else
                    {
                        countBad++;
                    }
                }
                if (countGood / (float)(countGood + countBad) < 0.6f) continue;

                informativeSegmentCount++;
                // Note the predicted purity for any ploidies that don't predict MAF=0.5:
                List<double> possiblePurities = new List<double>();
                for (int majorChromosomeCount = CN; majorChromosomeCount * 2 > CN; majorChromosomeCount--)
                {
                    double pureMAF = 1 - (majorChromosomeCount / (double)CN);
                    if (MedianMAF < pureMAF - 0.05f)
                    {
                        // This doesn't look plausible.  At 100% purity, we expect MedianMAF to be approximately equal
                        // to pureMAF (and MedianMAF rises to 0.5 as purity drops to 0%).  
                        continue;
                    }

                    double purity = (MedianMAF - 0.5) / (pureMAF - 0.5);
                    if (MedianMAF < pureMAF) purity = 1; // sanity-check
                    if (purity < 0.25) continue; // For now, *assume* purity is reasonably high
                    possiblePurities.Add(purity);
                }
                if (possiblePurities.Count > 1)
                {
                    continue;
                }
                foreach (double purity in possiblePurities)
                {
                    allPurities.Add(purity);
                    totalPurity += purity / possiblePurities.Count;
                }
                totalWeight++;
            }
            Console.WriteLine("Reviewed {0} informative segments (attempted segments {1})", informativeSegmentCount, attemptedSegmentCount);
            if (totalWeight == 0 || allPurities.Count == 0)
            {
                Console.Error.WriteLine("Warning: DerivePurityEstimateFromVF unable to model tumor purity from this data-set and truth set");
            }
            else
            {
                Console.WriteLine("Mean purity: {0:F4}", totalPurity / totalWeight);
                allPurities.Sort();
                Console.WriteLine("Median purity: {0:F4}", allPurities[allPurities.Count / 2]);
            }
            Console.WriteLine("DerivePurityEstimateFromVF complete.");
        }


        /// <summary>
        /// Compute the expected bin counts for each copy number, given a specified bin count for CN=2 regions, and assuming
        /// pure tumor data.  
        /// </summary>
        protected static double[] GetProjectedMeanCoverage(double diploidCoverage, int maximumCopyNumber)
        {
            double[] mu = new double[maximumCopyNumber + 1];
            for (int count = 0; count < mu.Length; count++)
            {
                mu[count] = diploidCoverage * count / 2f;
            }
            return mu;
        }

        /// <summary>
        ///  Initialize model points by subsampling from existing segment Coverage and MAF values. 
        ///  Use distanceThreshold to ensure that both large and small cluster components get subsampled
        /// </summary>
        protected List<ModelPoint> InitializeModelPoints(List<SegmentInfo> segments, int numClusters, double distanceThreshold)
        {
            List<ModelPoint> modelPoints = new List<ModelPoint>();
            List<SegmentInfo> usableSegments = new List<SegmentInfo>();
            List<SegmentInfo> usedSegments = new List<SegmentInfo>();

            foreach (SegmentInfo segment in segments)
            {
                if (segment.Cluster != -1 && segment.MAF >= 0)
                    usableSegments.Add(segment);
            }

            Random rnd = new Random();
            int lastIndex = rnd.Next(1, usableSegments.Count);
            usedSegments.Add(usableSegments[lastIndex]);
            int counter = 1;
            double attempts = 0;
            while (counter < numClusters)
            {
                int newIndex = rnd.Next(1, usableSegments.Count);
                attempts += 1.0;
                double distance = GetModelDistance(usableSegments[lastIndex].Coverage, usableSegments[newIndex].Coverage, usableSegments[lastIndex].MAF, usableSegments[newIndex].MAF);
                if (distance > distanceThreshold || attempts / usableSegments.Count > 0.3) // escape outlier minima
                {
                    usedSegments.Add(usableSegments[newIndex]);
                    counter++;
                    lastIndex = newIndex;
                    attempts = 0;
                }
            }
            // Initialize model points with coverage and MAF values from subsampled segments
            for (int i = 0; i < numClusters; i++)
            {
                ModelPoint point = new ModelPoint();
                point.Coverage = usedSegments[i].Coverage;
                point.MAF = usedSegments[i].MAF;
                SegmentPloidy ploidy = new SegmentPloidy();
                ploidy.CopyNumber = 2;
                ploidy.MajorChromosomeCount = 1;
                point.Ploidy = ploidy;
                point.Cluster = i + 1;
                modelPoints.Add(point);
            }
            return modelPoints;
        }

        /// <summary>
        ///  Initialize model points given diploid purity modelInitialize model points given somatic purity model
        /// </summary>
        protected List<ModelPoint> InitializeModelPoints(List<SegmentInfo> segments, double coverage, int percentPurity, int numClusters)
        {
            List<ModelPoint> modelPoints = new List<ModelPoint>();
            SomaticCaller.CoveragePurityModel model = new SomaticCaller.CoveragePurityModel(somaticCallerParameters.MaximumCopyNumber);
            model.DiploidCoverage = coverage;
            model.Purity = percentPurity / 100f;

            double[] mu = GetProjectedMeanCoverage(model.DiploidCoverage, somaticCallerParameters.MaximumCopyNumber);
            double diploidMAF = this.AllPloidies[3].MinorAlleleFrequency; /// %%% Magic number!

            /////////////////////////////////////////////
            // Update the parameters in each SegmentPloidy object, and construct corresponding SegmentInfo objects
            foreach (SegmentPloidy ploidy in this.AllPloidies)
            {
                ModelPoint point = new ModelPoint();
                double pureCoverage = mu[ploidy.CopyNumber];
                point.Coverage = (model.Purity * pureCoverage) + (1 - model.Purity) * model.DiploidCoverage;
                double pureMAF = ploidy.MinorAlleleFrequency;
                if (ploidy.MajorChromosomeCount * 2 == ploidy.CopyNumber)
                {
                    point.MAF = (model.Purity * ploidy.CopyNumber * pureMAF) + ((1 - model.Purity) * 2 * diploidMAF);
                    point.MAF /= model.Purity * ploidy.CopyNumber + (1 - model.Purity) * 2;
                    if (double.IsNaN(point.MAF)) point.MAF = 0;
                }
                else
                {
                    point.MAF = (model.Purity * ploidy.CopyNumber * pureMAF) + ((1 - model.Purity) * 1);
                    point.MAF /= model.Purity * ploidy.CopyNumber + (1 - model.Purity) * 2;
                }
                point.Ploidy = ploidy;
                modelPoints.Add(point);
                point.CN = ploidy.CopyNumber;
                ploidy.MixedMinorAlleleFrequency = point.MAF;
                ploidy.MixedCoverage = point.Coverage;
            }

            // estimate distance between each model point and segments 
            List<double> modelPointsScore = new List<double>();
            foreach (ModelPoint modelPoint in modelPoints)
            {
                List<double> distanceList = new List<double>();
                foreach (SegmentInfo info in segments)
                {
                    if (info.MAF >= 0)
                        distanceList.Add(GetModelDistance(info.Coverage, modelPoint.Coverage, info.MAF, modelPoint.MAF));
                }
                distanceList.Sort();
                double v15th_percentile = distanceList[Convert.ToInt32(distanceList.Count * 0.15)];
                // use model points with good fit to observed values
                modelPointsScore.Add(v15th_percentile);
            }
            // sort list and return indices
            var sortedScores = modelPointsScore.Select((x, i) => new KeyValuePair<double, int>(x, i)).OrderBy(x => x.Key).ToList();
            List<double> scoresValue = sortedScores.Select(x => x.Key).ToList();
            List<int> scoresIndex = sortedScores.Select(x => x.Value).ToList();

            List<ModelPoint> selectedModelPoints = new List<ModelPoint>();

            for (int i = 0; i < numClusters; i++)
            {
                modelPoints[scoresIndex[i]].Cluster = i + 1;
                selectedModelPoints.Add(modelPoints[scoresIndex[i]]);
            }

            return selectedModelPoints;
        }


        // Initialize model points given expected ploidy and purity values 
        protected List<ModelPoint> InitializeModelPoints(SomaticCaller.CoveragePurityModel model)
        {
            List<ModelPoint> modelPoints = new List<ModelPoint>();

            double[] mu = GetProjectedMeanCoverage(model.DiploidCoverage, somaticCallerParameters.MaximumCopyNumber);
            double diploidMAF = this.AllPloidies[3].MinorAlleleFrequency; /// %%% Magic number!

            /////////////////////////////////////////////
            // Update the parameters in each SegmentPloidy object, and construct corresponding SegmentInfo objects
            foreach (SegmentPloidy ploidy in this.AllPloidies)
            {
                ModelPoint point = new ModelPoint();
                double pureCoverage = mu[ploidy.CopyNumber];
                point.Coverage = (model.Purity * pureCoverage) + (1 - model.Purity) * model.DiploidCoverage;
                double pureMAF = ploidy.MinorAlleleFrequency;
                if (ploidy.MajorChromosomeCount * 2 == ploidy.CopyNumber)
                {
                    point.MAF = (model.Purity * ploidy.CopyNumber * pureMAF) + ((1 - model.Purity) * 2 * diploidMAF);
                    point.MAF /= model.Purity * ploidy.CopyNumber + (1 - model.Purity) * 2;
                    if (double.IsNaN(point.MAF)) point.MAF = 0;
                }
                else
                {
                    point.MAF = (model.Purity * ploidy.CopyNumber * pureMAF) + ((1 - model.Purity) * 1);
                    point.MAF /= model.Purity * ploidy.CopyNumber + (1 - model.Purity) * 2;
                }
                point.Ploidy = ploidy;
                modelPoints.Add(point);
                point.CN = ploidy.CopyNumber;
                ploidy.MixedMinorAlleleFrequency = point.MAF;
                ploidy.MixedCoverage = point.Coverage;
            }

            return modelPoints;
        }

        /// <summary>
        /// Fit a Gaussian mixture model.
        /// Fix the means to the model MAF and Coverage and run the EM algorithm until convergence.
        /// Compute the empirical MAF and Coverage.
        /// Fix the means to the empirical MAF and Coverage and run the EM algorithm again until convergence.
        /// Always estimate the full covariance matrix?
        /// </summary>
        /// <param name="model"></param>
        /// <param name="segments"></param>
        /// <param name="debugPath"></param>
        /// <returns></returns>
        protected double FitGaussians(SomaticCaller.CoveragePurityModel model, List<SegmentInfo> segments, string debugPath = null, double knearestNeighbourCutoff = Int32.MaxValue)
        {
            List<ModelPoint> modelPoints = InitializeModelPoints(model);

            GaussianMixtureModel gmm = new GaussianMixtureModel(modelPoints, segments, this.MeanCoverage, this.CoverageWeightingFactor, knearestNeighbourCutoff);
            double likelihood = gmm.Fit();

            if (debugPath != null)
            {
                // write Gaussian mixture model to debugPath
                using (StreamWriter writer = new StreamWriter(debugPath))
                {
                    writer.WriteLine("CN\tMajor Chr #\tMAF\tCoverage\tOmega\tMu0\tMu1\tSigma00\tSigma01\tSigma10\tSigma11");
                    foreach (ModelPoint modelPoint in modelPoints)
                    {
                        writer.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}",
                            modelPoint.Ploidy.CopyNumber, modelPoint.Ploidy.MajorChromosomeCount,
                            modelPoint.Ploidy.MixedMinorAlleleFrequency, modelPoint.Ploidy.MixedCoverage,
                            modelPoint.Ploidy.Omega, modelPoint.Ploidy.Mu[0], modelPoint.Ploidy.Mu[1],
                            modelPoint.Ploidy.Sigma[0][0], modelPoint.Ploidy.Sigma[0][1],
                            modelPoint.Ploidy.Sigma[1][0], modelPoint.Ploidy.Sigma[1][1]);
                    }

                    writer.WriteLine("");
                    writer.WriteLine("MAF\tCoverage\tPosterior Probabilities");
                    StringBuilder sb = new StringBuilder();
                    foreach (SegmentInfo segment in segments)
                    {
                        sb.Clear();
                        sb.AppendFormat("{0}\t{1}", segment.MAF, segment.Coverage);
                        foreach (ModelPoint modelPoint in modelPoints)
                        {
                            sb.AppendFormat("\t{0}", segment.PosteriorProbs[modelPoint]);
                        }
                        writer.WriteLine(sb.ToString());
                    }
                }
            }

            return likelihood;
        }

        /// <summary>
        /// Given genome-wide copy number (CN) profile of the model estimate the total number of rearrangements that 
        /// need to be applied to a diploid genome to transform it into the tumor genome under given purity model. 
        /// The following logic is used:
        ///     1)	Assign one rearrangement score to a single CN state transition, i.e. transition 2 -> 3 will get a score of one 
        ///     while transition 2 -> 4 will get a score of 2. 
        ///     2)	Cumulative PercentCN of 80% and more for copy number bins > 2 indicate possible genome doubling. 
        ///     Assign score of 1 for genome doubling event. Use copy number 4 baseline instead of 2 and count events as in step 1.
        /// </summary>
        protected double DiploidModelDistance(SomaticCaller.CoveragePurityModel model, List<SegmentInfo> usableSegments, long genomeLength)
        {
            double totalCNevents = 0;
            int modelBaseline = 2;
            double amplificationPercentCN = 0;
            for (int copyNumber = 3; copyNumber < somaticCallerParameters.MaximumCopyNumber; copyNumber++)
                amplificationPercentCN += model.PercentCN[copyNumber];
            if (amplificationPercentCN > 0.8)
            {
                modelBaseline = 4;
                totalCNevents += 1;
            }
            for (int i = 0; i < model.CNs.Count; i++)
            {
                totalCNevents += Math.Abs(model.CNs[i] - modelBaseline) * (usableSegments[i].Segment.End - usableSegments[i].Segment.Begin) / (double)genomeLength;
            }
            model.DiploidDistance = (double)1.0 / Math.Max(0.001, totalCNevents);
            return totalCNevents;
        }

        /// <summary>
        /// Estimate genome distance between two purity models (weighted absolute difference between copy number profiles)
        /// /// </summary>
        protected double CalculateModelDistance(SomaticCaller.CoveragePurityModel model1, SomaticCaller.CoveragePurityModel model2, List<SegmentInfo> usableSegments, long genomeLength)
        {
            double genomeDistance = 0;
            // every model should have the same number of segments
            if (model1.CNs.Count != model2.CNs.Count)
            {
                Console.WriteLine("Models do not have the same number of usable CN segments");
                return 1;
            }
            for (int i = 0; i < model1.CNs.Count; i++)
            {
                genomeDistance += Math.Abs(model1.CNs[i] - model2.CNs[i]) * (usableSegments[i].Segment.End - usableSegments[i].Segment.Begin) / (double)genomeLength;
            }
            return genomeDistance;
        }

        /// <summary>
        /// Return the squared euclidean distance between (coverage, maf) and (coverage2, maf2) in scaled coverage/MAF space.
        /// </summary>
        protected double GetModelDistance(double coverage, double coverage2, double? maf, double maf2)
        {
            double diff = (coverage - coverage2) * CoverageWeightingFactor;
            double distance = diff * diff;
            if (!maf.HasValue || maf < 0) return 2 * distance;
            diff = (double)maf - maf2;
            distance += diff * diff;
            return distance;
        }

        /// <summary>
        /// Compute Silhouette coefficient https://en.wikipedia.org/wiki/Silhouette_(clustering)
        /// </summary>
        protected double ComputeSilhouette(List<SegmentInfo> usableSegments, int numClusters)
        {
            List<List<double>> withinClusterDistance = new List<List<double>>();
            List<List<double>> betweenClusterDistance = new List<List<double>>();
            for (int i = 0; i < numClusters; i++) withinClusterDistance.Add(new List<double>());
            for (int i = 0; i < numClusters; i++) betweenClusterDistance.Add(new List<double>());

            for (int k = 0; k < numClusters; k++)
            {
                for (int i = 0; i < usableSegments.Count; i++)
                {
                    for (int j = 0; j < usableSegments.Count; j++)
                    {
                        if (i != j && usableSegments[i].Cluster != -1 && usableSegments[j].Cluster != -1 && usableSegments[i].Cluster == k + 1 && usableSegments[i].MAF >= 0 && usableSegments[j].MAF >= 0)
                        {
                            if (usableSegments[i].Cluster == usableSegments[j].Cluster)
                            {
                                withinClusterDistance[k].Add(GetModelDistance(usableSegments[i].Coverage, usableSegments[j].Coverage, usableSegments[i].MAF, usableSegments[j].MAF));
                            }
                            else
                            {
                                betweenClusterDistance[k].Add(GetModelDistance(usableSegments[i].Coverage, usableSegments[j].Coverage, usableSegments[i].MAF, usableSegments[j].MAF));
                            }
                        }
                    }
                }
            }
            double silhouette = 0;
            for (int i = 0; i < numClusters; i++)
            {
                if (withinClusterDistance[i].Count > 2 && betweenClusterDistance[i].Count > 2)
                {
                    double a = CanvasCommon.Utilities.Median(withinClusterDistance[i]);
                    double b = CanvasCommon.Utilities.Median(betweenClusterDistance[i]);
                    silhouette += (b - a) / Math.Max(a, b);
                }
            }
            return silhouette / numClusters;
        }

        /// <summary>
        /// Refine our estimate of diploid MAF.  Our baseline model, based on a binomial distribution, doesn't match real-world data 
        /// perfectly - and it seems some real-world samples have higher or lower diploid MAF than others even for the same coverage level.
        /// So, we average together the initial model with an empirical fit.
        /// </summary>
        protected void RefineDiploidMAF(List<SegmentInfo> segments, List<ModelPoint> modelPoints)
        {
            // First pass: Find segments assigned to even copy number (no LOH), get the mean MAF, and use that to refine our MAF model:
            double[] diploidMAF = new double[1 + somaticCallerParameters.MaximumCopyNumber / 2];
            double[] diploidMAFWeight = new double[1 + somaticCallerParameters.MaximumCopyNumber / 2];

            // Seed the model with 10 megabases of coverage at the "expected" MAF, so that we don't end up with silly values
            // based on one or two point:
            double dummyWeight = 10000000;
            foreach (ModelPoint modelPoint in modelPoints)
            {
                if (modelPoint.CN % 2 == 1) continue;
                if (modelPoint.Ploidy.MajorChromosomeCount * 2 != modelPoint.CN) continue;
                diploidMAF[modelPoint.CN / 2] += dummyWeight * modelPoint.MAF;
                diploidMAFWeight[modelPoint.CN / 2] += dummyWeight;
            }

            foreach (SegmentInfo info in segments)
            {
                if (info.MAF < 0) continue;
                double bestDistance = double.MaxValue;
                int bestCN = 0;
                ModelPoint bestModelPoint = null;
                foreach (ModelPoint modelPoint in modelPoints)
                {
                    double distance = GetModelDistance(info.Coverage, modelPoint.Coverage, info.MAF, modelPoint.MAF);

                    if (distance < bestDistance)
                    {
                        bestDistance = distance;
                        bestCN = modelPoint.CN;
                        info.Ploidy = modelPoint.Ploidy;
                        bestModelPoint = modelPoint;
                    }
                }
                if (bestModelPoint.CN % 2 == 0 && bestModelPoint.Ploidy.MajorChromosomeCount * 2 == bestModelPoint.CN)
                {
                    if (info.MAF < 0.4) continue;
                    diploidMAF[bestModelPoint.CN / 2] += info.Weight * info.MAF;
                    diploidMAFWeight[bestModelPoint.CN / 2] += info.Weight;
                }
            }

            foreach (ModelPoint modelPoint in modelPoints)
            {
                if (modelPoint.CN % 2 == 1) continue;
                if (modelPoint.Ploidy.MajorChromosomeCount * 2 != modelPoint.CN) continue;
                modelPoint.MAF = diploidMAF[modelPoint.CN / 2] / diploidMAFWeight[modelPoint.CN / 2];
            }
        }


        /// <summary>
        /// Helper function for ModelDeviation. Outputs
        /// estimates of average cluster deviation.
        /// </summary>
        protected double ClusterDeviation(SomaticCaller.CoveragePurityModel model, List<ModelPoint> modelPoints, List<double> centroidMAFs, List<double> centroidCoverage, List<SegmentInfo> segments, int numClusters, double tempDeviation, out int heterogeneousClusters, out double heterogeneityIndex, bool bestModel, string debugPathClusterInfo = null)
        {
            // compute average deviation for each cluster (clusterDeviation)
            List<ClusterModel> clusterDeviations = new List<ClusterModel>();
            for (int clusterID = 0; clusterID < numClusters; clusterID++)
            {
                ClusterModel clusterInfo = new ClusterModel();
                clusterInfo.ClusterID = clusterID + 1;
                clusterDeviations.Add(clusterInfo);
            }

            // populate clusters
            List<List<double>> clusterDistance = new List<List<double>>();
            for (int i = 0; i < numClusters; i++) clusterDistance.Add(new List<double>());
            List<double> tmpClusterDistance = new List<double>();

            foreach (ModelPoint modelPoint in modelPoints)
            {
                if (modelPoint.Coverage < MeanCoverage * 2.0)
                {
                    double bestDeviation = Double.MaxValue;
                    double currentDeviation = Double.MaxValue;
                    int bestCluster = 0;
                    for (int i = 0; i < numClusters; i++)
                    {
                        currentDeviation = GetModelDistance(centroidCoverage[i], modelPoint.Coverage, centroidMAFs[i], modelPoint.MAF);
                        if (currentDeviation < bestDeviation)
                        {
                            bestCluster = i;
                            bestDeviation = currentDeviation;
                        }
                    }
                    clusterDistance[bestCluster].Add(bestDeviation);
                    tmpClusterDistance.Add(bestDeviation);
                }
            }

            double clusterDeviation = 0;
            double clusterOverenrichment = 0;
            for (int i = 0; i < numClusters; i++)
            {
                clusterDeviation += clusterDistance[i].Sum();
                clusterOverenrichment += clusterDistance[i].Count();
            }
            clusterDeviation /= numClusters;
            heterogeneityIndex = clusterOverenrichment / numClusters;

            if (tmpClusterDistance.Count == 0)
            {
                heterogeneousClusters = int.MaxValue;
                return Double.MaxValue;
            }

            double medianClusterDistance = CanvasCommon.Utilities.Median(tmpClusterDistance);

            // compute cluster mean and standard deviation
            for (int clusterID = 0; clusterID < numClusters; clusterID++)
            {
                if (clusterDistance[clusterID].Count > 1)
                    clusterDeviations[clusterID].ClusterMedianDistance = CanvasCommon.Utilities.Median(clusterDistance[clusterID]);
                else if (clusterDistance[clusterID].Count == 1)
                    clusterDeviations[clusterID].ClusterMedianDistance = clusterDistance[clusterID][0];
                else
                    clusterDeviations[clusterID].ClusterMedianDistance = tmpClusterDistance.Max();
            }

            // exlcude clusters with deviation larger than 1.25 of average deviation
            // these clusters locate far from expected model centroids and most likely represent segments coming from heterogeneous variants 
            List<double> heterogeneousClusterID = new List<double>();
            foreach (ClusterModel clusterInfo in clusterDeviations)
                if (clusterInfo.ClusterMedianDistance > medianClusterDistance * somaticCallerParameters.HeterogeneousClusterMedianCutoff)
                    heterogeneousClusterID.Add((int)clusterInfo.ClusterID);


            // store signatures of potential heterogeneous variants 
            if (heterogeneousClusterID.Count > 0 && bestModel)
            {
                foreach (SegmentInfo info in segments)
                    if (heterogeneousClusterID.Contains(info.Cluster.Value))
                        this.HeterogeneousSegmentsSignature.Add(info.Segment.Begin + info.Segment.End + info.Segment.Counts.Count);
            }
            heterogeneousClusters = heterogeneousClusterID.Count;


            //  write cluster deviations
            if (debugPathClusterInfo != null)
            {
                using (StreamWriter debugWriter = new StreamWriter(debugPathClusterInfo))
                {
                    // Write clustering results
                    debugWriter.WriteLine("#MAF\tCoverage\tClusterID\tClusterMedianDistance\tClusterDistanceIQR");
                    for (int i = 0; i < segments.Count; i++)
                    {
                        debugWriter.Write("{0}\t{1}\t{2}\t", segments[i].MAF, segments[i].Coverage, segments[i].Cluster);
                        if (segments[i].Cluster.HasValue && segments[i].Cluster > 0)
                        {
                            int cluster = (int)segments[i].Cluster - 1;
                            debugWriter.Write("{0}", clusterDeviations[cluster].ClusterMedianDistance);
                        }
                        debugWriter.WriteLine();
                    }
                }
            }
            return clusterDeviation;
        }

        /// <summary>
        /// Helper function for ModelOverallCoverageAndPurity.  Measure the deviation (mismatch) between our
        /// model of expected coverage + minor allele frequency, and the actual data.
        /// Note that this method updates the parameters in this.AllPloidies to match this model.
        /// TotalDeviation = PrecisionWeight * PrecisionDeviation + (1 - PrecisionWeight) * AccuracyDeviation
        /// PrecisionDeviation is the weighted average of the distance between segments and their assigned ploidy
        /// AccuracyDeviation is the weighted average of the distance from the segment centroid 
        /// and the corresponding ploidy.
        /// </summary>
        protected double ModelDeviation(List<double> centroidMAFs, List<double> centroidCoverage, SomaticCaller.CoveragePurityModel model, List<SegmentInfo> segments, int numClusters, string debugPathClusterInfo = null, bool bestModel = false, string debugPath = null)
        {
            List<ModelPoint> modelPoints = InitializeModelPoints(model);
            double precisionDeviation = 0;
            this.RefineDiploidMAF(segments, modelPoints);

            /////////////////////////////////////////////
            // Cluster our segments:
            Array.Clear(model.PercentCN, 0, model.PercentCN.Length);
            model.CNs.Clear();
            double totalWeight = 0;
            double totalBasesNormal = 0;
            foreach (SegmentInfo info in segments)
            {
                double bestDistance = double.MaxValue;
                int bestCN = 0;
                ModelPoint bestModelPoint = null;
                foreach (ModelPoint modelPoint in modelPoints)
                {
                    double distance = GetModelDistance(info.Coverage, modelPoint.Coverage, info.MAF, modelPoint.MAF);
                    if (distance < bestDistance)
                    {
                        bestDistance = distance;
                        bestCN = modelPoint.CN;
                        info.Ploidy = modelPoint.Ploidy;
                        bestModelPoint = modelPoint;
                    }
                }

                bestDistance = Math.Sqrt(bestDistance);
                info.Distance = bestDistance;
                precisionDeviation += bestDistance * info.Weight;
                totalWeight += info.Weight;
                model.PercentCN[bestCN] += info.Weight;
                if (bestCN == 2 && info.Ploidy.MajorChromosomeCount == 1) totalBasesNormal += info.Weight;
                bestModelPoint.Weight += info.Weight;
                bestModelPoint.EmpiricalCoverage += info.Weight * info.Coverage;
                if (info.MAF >= 0)
                {
                    bestModelPoint.EmpiricalMAF += info.Weight * info.MAF;
                    bestModelPoint.MAFWeight += info.Weight;
                }
                // add CN variant of the segment to the model 
                if (bestCN == 2 && info.Ploidy.MajorChromosomeCount == 2)
                    // aproximate LOH; we presume that LOH counts as one event, hence similar in effect to HET deletions
                    model.CNs.Add(1);
                else
                    model.CNs.Add(bestCN);
            }
            precisionDeviation /= totalWeight;

            // Compute AccuracyDeviation:
            double accuracyDeviation = 0;
            foreach (ModelPoint modelPoint in modelPoints)
            {
                if (modelPoint.Weight == 0) continue;
                modelPoint.EmpiricalCoverage /= modelPoint.Weight;
                if (modelPoint.MAFWeight > 0) modelPoint.EmpiricalMAF /= modelPoint.MAFWeight;
                double distance = this.GetModelDistance(modelPoint.Coverage, modelPoint.EmpiricalCoverage, modelPoint.MAF, modelPoint.EmpiricalMAF);
                distance = Math.Sqrt(distance);
                accuracyDeviation += distance * modelPoint.Weight;
                if (!string.IsNullOrEmpty(debugPath))
                {
                    Console.WriteLine("{0}\t{1}\t{2:F2}\t{3:F0}\t{4:F2}\t{5:F0}\t{6:F3},{7:F0}",
                        modelPoint.CN, modelPoint.Ploidy.MajorChromosomeCount,
                        modelPoint.MAF, modelPoint.Coverage,
                        modelPoint.EmpiricalMAF, modelPoint.EmpiricalCoverage,
                        distance, modelPoint.Weight);
                }
            }
            accuracyDeviation /= totalWeight;

            // estimate abundance of each CN state
            for (int index = 0; index < model.PercentCN.Length; index++)
            {
                model.PercentCN[index] /= totalWeight;
            }

            // get model ploidy
            for (int index = 0; index < model.PercentCN.Length; index++)
            {
                model.Ploidy += index * model.PercentCN[index];
            }

            // standard somatic model deviation
            double tempDeviation = precisionDeviation * 0.5f + 0.5f * accuracyDeviation;

            // enough segments to compute cluster deviation?
            int heterogeneousClusters = 0;
            double heterogeneityIndex = 0;
            double clusterDeviation = 0;
            int validMAFCount = segments.Count(x => x.MAF >= 0);
            if (validMAFCount > 100 && segments.Count > 100 && centroidMAFs.Count < 10 && !this.IsEnrichment)
            {
                // compute cluster deviation
                clusterDeviation = ClusterDeviation(model, modelPoints, centroidMAFs, centroidCoverage, segments, numClusters, tempDeviation, out heterogeneousClusters, out heterogeneityIndex, bestModel, debugPathClusterInfo);
            }

            // compute total deviation
            double totalDeviation;
            if (heterogeneousClusters > somaticCallerParameters.HeterogeneousClustersCutoff)
                totalDeviation = somaticCallerParameters.PrecisionWeightingFactor * precisionDeviation + somaticCallerParameters.PrecisionWeightingFactor * accuracyDeviation + somaticCallerParameters.PrecisionWeightingFactor * clusterDeviation;
            else
                totalDeviation = tempDeviation;


            model.PercentNormal = totalBasesNormal / totalWeight;
            if (!string.IsNullOrEmpty(debugPath))
            {
                try
                {
                    using (StreamWriter debugWriter = new StreamWriter(debugPath))
                    {
                        debugWriter.WriteLine("#MAF\tCoverage\t");
                        foreach (ModelPoint modelPoint in modelPoints)
                        {
                            string gt = modelPoint.Ploidy.MajorChromosomeCount.ToString() + "/" + modelPoint.CN.ToString();
                            debugWriter.WriteLine("{0}\t{1}\t{2}\t", modelPoint.MAF, modelPoint.Coverage, gt);
                        }
                        debugWriter.WriteLine();
                        debugWriter.WriteLine("#MAF\tCoverage\tBestDistance\tChromosome\tBegin\tEnd\tLength\tTruthSetCN");
                        foreach (SegmentInfo info in segments)
                        {
                            // Find the best fit for this segment:
                            double bestDistance = double.MaxValue;
                            foreach (ModelPoint modelPoint in modelPoints)
                            {
                                double distance = GetModelDistance(info.Coverage, modelPoint.Coverage, info.MAF, modelPoint.MAF);
                                if (distance < bestDistance) bestDistance = distance;
                            }
                            bestDistance = Math.Sqrt(bestDistance);
                            debugWriter.Write("{0}\t{1}\t", info.MAF, info.Coverage);
                            debugWriter.Write("{0}\t{1}\t{2}\t{3}\t", bestDistance, info.Segment.Chr, info.Segment.Begin, info.Segment.End);
                            debugWriter.Write("{0}\t", info.Segment.End - info.Segment.Begin);
                            int CN = this.GetKnownCNForSegment(info.Segment);
                            debugWriter.Write("{0}\t", CN);
                            debugWriter.WriteLine();
                        }
                    }
                }
                catch (IOException ex)
                {
                    // Whine, but continue - not outputing this file is not fatal.
                    Console.Error.WriteLine(ex.ToString());
                }
            }

            // make sure that CN profile length is equal to the usable segments length
            if (model.CNs.Count != segments.Count)
            {
                throw new IndexOutOfRangeException(String.Concat("Canvas Somatic Caller error: index sizes do not match, ",
                    model.CNs.Count, " != ", segments.Count));
            }
            model.PrecisionDeviation = precisionDeviation;
            model.AccuracyDeviation = accuracyDeviation;
            model.Deviation = totalDeviation;
            model.HeterogeneityIndex = heterogeneityIndex;
            model.ClusterDeviation = clusterDeviation;
            return totalDeviation;
        }

        protected class CoveragePurityModel : CoverageModel
        {
            public double Purity;
            public double PercentNormal; // Percentage of bases called as CN=2, MajorChromosomeCount=1 (no LOH)
            public double DiploidDistance;
            public double? InterModelDistance;
            public double? HeterogeneityIndex;
            public double? ClusterDeviation;
            public double[] PercentCN;

            public CoveragePurityModel(int maximumCopyNumber)
            {
                PercentCN = new double[maximumCopyNumber + 1];
            }


        }

        static public List<SegmentInfo> GetUsableSegmentsForModeling(List<CanvasSegment> segments, bool IsEnrichment, int minimumVariantFrequenciesForInformativeSegment)
        {
            // Get the average count everwhere.  Exclude segments whose coverage is >2x this average.
            List<float> tempCountsList = new List<float>();
            foreach (CanvasSegment segment in segments)
            {
                if (IsEnrichment)
                {
                    tempCountsList.Add(Convert.ToSingle(CanvasCommon.Utilities.Median(segment.Counts)));
                }
                else
                {
                    foreach (float value in segment.Counts)
                    {
                        tempCountsList.Add(value);
                    }
                }
            }

            Tuple<float, float, float> coverageQuartiles = CanvasCommon.Utilities.Quartiles(tempCountsList);
            float overallMedian = coverageQuartiles.Item2;

            List<SegmentInfo> usableSegments = new List<SegmentInfo>();
            foreach (CanvasSegment segment in segments)
            {
                if (segment.End - segment.Begin < 5000) continue;
                SegmentInfo info = new SegmentInfo();
                info.Segment = segment;
                // If the segment has few or no variants, then don't use the MAF for this segment - set to -1 (no frequency)
                // Typically a segment will have no variants if it's on chrX or chrY and starling knows not to call a
                // heterozygous variant there (other than in the PAR regions).
                if (segment.VariantFrequencies.Count < minimumVariantFrequenciesForInformativeSegment)
                {
                    info.MAF = -1;
                }
                else
                {
                    List<double> MAF = new List<double>();
                    foreach (float value in segment.VariantFrequencies) MAF.Add(value > 0.5 ? 1 - value : value);
                    MAF.Sort();
                    info.MAF = MAF[MAF.Count / 2];
                }
                info.Coverage = CanvasCommon.Utilities.Median(segment.Counts);
                if (info.Coverage > overallMedian * 2) continue;
                if (segments.Count > 100)
                {
                    info.Weight = segment.End - segment.Begin;
                }
                else
                {
                    info.Weight = segment.BinCount;
                }
                if (segment.VariantFrequencies.Count < 10)
                {
                    info.Weight *= (double)segment.VariantFrequencies.Count / 10;
                }
                usableSegments.Add(info);
            }
            return usableSegments;
        }


        /// <summary>
        /// Estimate the optimal number of clusters in an Expectation Maximizationfrom model using silhouette coefficient
        /// </summary>
        public List<ModelPoint> BestNumClusters(List<SegmentInfo> usableSegments, double medianCoverageLevel, double bestCoverageWeightingFactor, double knearestNeighbourCutoff)
        {
            double bestSilhouette = double.MinValue;
            List<ModelPoint> bestModelPoints = new List<ModelPoint>();
            int bestNumClusters = 0;
            int maxNumClusters = 8;

            // find distanceThreshold, use it in InitializeModelPoints
            List<double> tempModelDistanceList = new List<double>();
            for (int i = 0; i < usableSegments.Count; i++)
            {
                for (int j = 0; j < usableSegments.Count; j++)
                {
                    if (i != j && usableSegments[i].Cluster != -1 && usableSegments[j].Cluster != -1 && usableSegments[i].MAF >= 0 && usableSegments[j].MAF >= 0)
                    {
                        tempModelDistanceList.Add(GetModelDistance(usableSegments[i].Coverage, usableSegments[j].Coverage, usableSegments[i].MAF, usableSegments[j].MAF));
                    }
                }
            }
            tempModelDistanceList.Sort();
            int distanceThresholdIndex = Math.Min(Convert.ToInt32(tempModelDistanceList.Count * 0.8), tempModelDistanceList.Count - 1);
            double distanceThreshold = tempModelDistanceList[distanceThresholdIndex]; // enable capturing clusters with less than 30% abundance 

            // find optimal number of clusters using silhouette distance and return bestModelPoints
            for (int numClusters = 4; numClusters < maxNumClusters; numClusters++)
            {
                for (int i = 0; i < 10; i++)
                {
                    List<ModelPoint> tempModelPoints = InitializeModelPoints(usableSegments, numClusters, distanceThreshold);
                    GaussianMixtureModel tempgmm = new GaussianMixtureModel(tempModelPoints, usableSegments, medianCoverageLevel, bestCoverageWeightingFactor, knearestNeighbourCutoff);
                    double currentLikelihood = tempgmm.runExpectationMaximization(); // return BIC rather than raw likelihood
                    double currentSilhouette = ComputeSilhouette(usableSegments, numClusters);
                    if (bestSilhouette < currentSilhouette)
                    {
                        bestSilhouette = currentSilhouette;
                        bestNumClusters = numClusters;
                        if (bestModelPoints.Count > 0)
                            bestModelPoints.Clear();
                        foreach (ModelPoint tempModelPoint in tempModelPoints)
                            bestModelPoints.Add(tempModelPoint);
                    }
                }
            }
            return bestModelPoints;
        }

        /// <summary>
        /// Find K-nearest neighbours cutoff. The cutoff is used to identify outlier segments that do not belong to any cluster
        /// </summary>
        public double KnearestNeighbourCutoff(List<SegmentInfo> usableSegments)
        {
            int kneighours = 10;
            List<double> knearestNeighbourList = new List<double>();
            for (int i = 0; i < usableSegments.Count; i++)
            {
                List<double> tempModelDistanceList = new List<double>();

                for (int j = 0; j < usableSegments.Count; j++)
                {
                    if (i != j)
                    {
                        tempModelDistanceList.Add(GetModelDistance(usableSegments[i].Coverage, usableSegments[j].Coverage, usableSegments[i].MAF, usableSegments[j].MAF));
                    }
                }
                tempModelDistanceList.Sort();
                double distance = 0;
                for (int k = 0; k < kneighours; k++)
                {
                    distance += tempModelDistanceList[k];
                }
                usableSegments[i].KnearestNeighbour = distance;
                knearestNeighbourList.Add(distance);
            }

            knearestNeighbourList.Sort();
            double knearestNeighbourCutoff = knearestNeighbourList[Convert.ToInt32(knearestNeighbourList.Count * 0.99)];
            return knearestNeighbourCutoff;
        }

        /// <summary>
        /// Estimate the optimal number of clusters in an Expectation Maximizationfrom model using silhouette coefficient
        /// </summary>
        public double BestCoverageWeightingFactor(List<SegmentInfo> usableSegments, int maxCoverageLevel, int medianCoverageLevel, double knearestNeighbourCutoff)
        {
            double bestLikelihood = double.MinValue;
            double bestCoverageWeightingFactor = 0;
            // magic-scaling for now - keep small to penalize coverage, tested on 50+ groundtruth corpus
            double maxCoverageWeightingFactor = somaticCallerParameters.CoverageWeighting / medianCoverageLevel;
            double minCoverageWeightingFactor = 0.1 / maxCoverageLevel;
            double stepCoverageWeightingFactor = Math.Max(0.00001, (maxCoverageWeightingFactor - minCoverageWeightingFactor) / 10);

            for (double coverageWeighting = minCoverageWeightingFactor; coverageWeighting < maxCoverageWeightingFactor; coverageWeighting += stepCoverageWeightingFactor)
            {
                List<ModelPoint> tempModelPoints = InitializeModelPoints(usableSegments, medianCoverageLevel / 2.0, 90, 6);
                GaussianMixtureModel tempgmm = new GaussianMixtureModel(tempModelPoints, usableSegments, medianCoverageLevel, coverageWeighting, knearestNeighbourCutoff);
                double currentLikelihood = tempgmm.runExpectationMaximization();
                if (currentLikelihood > bestLikelihood)
                {
                    bestLikelihood = currentLikelihood;
                    bestCoverageWeightingFactor = coverageWeighting;
                }
            }
            return bestCoverageWeightingFactor;
        }

        /// <summary>
        /// Identify the tuple (DiploidCoverage, OverallPurity) which best models our overall
        /// distribution of (MAF, Coverage) data across all segments.  Consider various tuples (first with a coarse-grained
        /// and then a fine-grained search), and for each one, measure the distortion - the average distance (weighted 
        /// by segment length) between actual and modeled (MAF, Coverage) coordinate.
        /// </summary>
        protected SomaticCaller.CoveragePurityModel ModelOverallCoverageAndPurity(long genomeLength, CanvasSomaticClusteringMode clusteringMode)
        {
            List<SegmentInfo> usableSegments;
            // Identify usable segments using our MinimumVariantFrequenciesForInformativeSegment cutoff, 
            // then (if we don't find enough) we can try again with progressively more permissive cutoffs.
            int validMAFCount = 0;
            while (true)
            {
                usableSegments = GetUsableSegmentsForModeling(this.Segments, IsEnrichment, somaticCallerParameters.MinimumVariantFrequenciesForInformativeSegment);
                validMAFCount = usableSegments.Count(x => x.MAF >= 0);
                if (validMAFCount > Math.Min(20, this.Segments.Count)) break; // We have enough usable segments with nonnull MAF
                if (somaticCallerParameters.MinimumVariantFrequenciesForInformativeSegment <= 5) break; // Give up on modeling
                somaticCallerParameters.MinimumVariantFrequenciesForInformativeSegment -= 15;
                somaticCallerParameters.MinimumVariantFrequenciesForInformativeSegment = Math.Max(5, somaticCallerParameters.MinimumVariantFrequenciesForInformativeSegment);
            }
            Console.WriteLine("Modeling overall coverage/purity across {0} segments", usableSegments.Count);
            if (usableSegments.Count < 10)
                throw new SomaticCaller.UncallableDataException("Cannot model coverage/purity with less than 10 segments.");

            // When computing distances between model and actual points, we want to provide roughly equal weight
            // to coverage (which covers a large range) and MAF, which falls in the range (0, 0.5).  
            // If we already knew the diploid coverage, then we'd know just how to scale things (catch-22).
            // Let's assume that the median coverage is a sane thing to use for scaling:
            List<float> tempCoverageList = new List<float>();
            List<double> knearestNeighbourList = new List<double>();

            // Segments clustering using Gaussian Expectation Maximisation

            // Step0: Prepare model parameters
            foreach (SegmentInfo info in usableSegments) tempCoverageList.Add(Convert.ToSingle(info.Coverage));
            Tuple<float, float, float> coverageQuartiles = CanvasCommon.Utilities.Quartiles(tempCoverageList);
            int minCoverageLevel = Convert.ToInt32(coverageQuartiles.Item1);
            int maxCoverageLevel = Convert.ToInt32(coverageQuartiles.Item3);
            int medianCoverageLevel = Convert.ToInt32(coverageQuartiles.Item2);
            this.CoverageWeightingFactor = somaticCallerParameters.CoverageWeighting / medianCoverageLevel;
            int bestNumClusters = 0;
            double knearestNeighbourCutoff = 0;
            List<double> centroidsMAF = new List<double>();
            List<double> centroidsCoverage = new List<double>();
            // Need  large number of segments for cluster analysis
            if (usableSegments.Count > 100 && validMAFCount > 100 && !this.IsEnrichment)
            {
                switch (clusteringMode)
                {
                    case CanvasSomaticClusteringMode.GaussianMixture:
                        // Step1: Find outliers
                        knearestNeighbourCutoff = KnearestNeighbourCutoff(usableSegments);

                        // Step2: Find the best CoverageWeightingFactor 
                        double bestCoverageWeightingFactor = BestCoverageWeightingFactor(usableSegments, maxCoverageLevel, medianCoverageLevel, knearestNeighbourCutoff);

                        // Step3: Find the optimal number of clusters
                        List<ModelPoint> modelPoints = BestNumClusters(usableSegments, medianCoverageLevel, bestCoverageWeightingFactor, knearestNeighbourCutoff);
                        bestNumClusters = modelPoints.Count;

                        // Step4: Find segment clusters using the final model
                        GaussianMixtureModel gmm = new GaussianMixtureModel(modelPoints, usableSegments, medianCoverageLevel, bestCoverageWeightingFactor, knearestNeighbourCutoff);
                        double likelihood = gmm.runExpectationMaximization();
                        break;

                    case CanvasSomaticClusteringMode.Density:
                        // Step1: Find outliers
                        knearestNeighbourCutoff = KnearestNeighbourCutoff(usableSegments);
                        // Step2: Density clustering 
                        double centroidCutoff = 0;
                        int clusterCount = 0;
                        List<int> numNumClusters = new List<int>();
                        List<double> centroidCutoffs = new List<double>();
                        double centoridStep = (somaticCallerParameters.UpperCentroidCutoff - somaticCallerParameters.LowerCentroidCutoff) / somaticCallerParameters.CentroidCutoffStep;
                        for (double centoridCutoff = somaticCallerParameters.LowerCentroidCutoff;
                            centoridCutoff < somaticCallerParameters.UpperCentroidCutoff;
                            centoridCutoff += centoridStep)
                            centroidCutoffs.Add(centoridCutoff);
                        centroidCutoffs.Reverse();

                        foreach (double centoridCutoff in centroidCutoffs)
                        {
                            DensityClusteringModel densityClustering = new DensityClusteringModel(usableSegments, CoverageWeightingFactor, knearestNeighbourCutoff, centoridCutoff);
                            densityClustering.EstimateDistance();
                            double distanceThreshold = densityClustering.EstimateDc();
                            densityClustering.GaussianLocalDensity(distanceThreshold);
                            densityClustering.FindCentroids();
                            clusterCount = densityClustering.FindClusters();
                            numNumClusters.Add(clusterCount);
                            Console.WriteLine(">>> Running density clustering for cutoff {0:F5} , number of clusters {1}", centoridCutoff, clusterCount);
                        }
                        var modeClustersValues = numNumClusters
                            .GroupBy(x => x)
                            .Select(g => new { Value = g.Key, Count = g.Count() })
                            .ToList(); // materialize the query to avoid evaluating it twice below
                        int maxCount = modeClustersValues.Max(g => g.Count); // throws InvalidOperationException if myArray is empty
                        IEnumerable<int> modes = modeClustersValues
                            .Where(g => g.Count == maxCount)
                            .Select(g => g.Value);
                        List<int> modesList = modes.ToList();
                        if (modesList.Count == 1)
                        {
                            clusterCount = modesList[0];
                            centroidCutoff = centroidCutoffs[numNumClusters.FindIndex(x => x == clusterCount)];

                        }
                        else if (modesList.Count == 2 || modesList.Count == 3)
                        {
                            if (modesList[1] < 7)
                                clusterCount = modesList[1];
                            else
                                clusterCount = modesList[0];
                            centroidCutoff = centroidCutoffs[numNumClusters.FindIndex(x => x == clusterCount)];
                        }
                        else
                        {
                            centroidCutoff = somaticCallerParameters.DefaultCentroidCutoff;
                            clusterCount = numNumClusters[centroidCutoffs.FindIndex(x => x == centroidCutoff)];

                        }

                        Console.WriteLine(">>> Running density selected cutoff {0:F5}", centroidCutoff);
                        DensityClusteringModel finalDc = new DensityClusteringModel(usableSegments, CoverageWeightingFactor, knearestNeighbourCutoff, centroidCutoff);
                        finalDc.EstimateDistance();
                        double finalDistanceThreshold = finalDc.EstimateDc();
                        finalDc.GaussianLocalDensity(finalDistanceThreshold);
                        finalDc.FindCentroids();
                        bestNumClusters = finalDc.FindClusters();
                        centroidsMAF = finalDc.GetCentroidsMAF();
                        centroidsCoverage = finalDc.GetCentroidsMAF();
                        break;
                    default:
                        throw new ApplicationException("Unsupported CanvasSomatic clustering mode: " + clusteringMode.ToString());
                }
            }

            // Note: Don't consider purity below 20 (at this point), becuase that creates a model that is very noise-sensitive.
            // We tried using a "tumor" sample that is actually just the real normal: We could overfit this data as very low 
            // purity 5% and make lots of (bogus) calls which fit the noise in coverage and MAF.

            double bestDeviation = double.MaxValue;
            List<SomaticCaller.CoveragePurityModel> allModels = new List<SomaticCaller.CoveragePurityModel>();
            // set best somatic model to pre-specified  ploidy and purity values
            if (this.userPloidy != null && this.userPurity != null)
            {
                SomaticCaller.CoveragePurityModel bestModel = new SomaticCaller.CoveragePurityModel(somaticCallerParameters.MaximumCopyNumber);
                bestModel.DiploidCoverage = GetDiploidCoverage(medianCoverageLevel, this.userPloidy.Value);
                bestModel.Purity = Convert.ToDouble(this.userPurity);

                this.ModelDeviation(centroidsMAF, centroidsCoverage, bestModel, usableSegments, bestNumClusters);
                this.DiploidModelDistance(bestModel, usableSegments, genomeLength);
                return bestModel;
            }
            // find best somatic model
            else
            {
                // Coarse search: Consider various (coverage, purity) tuples.  
                int minCoverage = (int)Math.Max(10, medianCoverageLevel / somaticCallerParameters.LowerCoverageLevelWeightingFactor);
                int maxCoverage = (int)Math.Max(10, medianCoverageLevel * somaticCallerParameters.UpperCoverageLevelWeightingFactor);
                int minPercentPurity = 20;
                int maxPercentPurity = 100;
                if (this.userPloidy != null)
                {
                    minCoverage = maxCoverage = (int)GetDiploidCoverage(medianCoverageLevel, this.userPloidy.Value);
                }
                if (this.userPurity != null)
                {
                    minPercentPurity = maxPercentPurity = (int)(this.userPurity.Value * 100);
                }
                int coverageStep = Math.Max(1, (maxCoverage - minCoverage) / somaticCallerParameters.CoverageLevelWeightingFactorLevels);
                Console.WriteLine(">>>DiploidCoverage: Consider {0}...{1} step {2}", minCoverage, maxCoverage, coverageStep);
                for (int coverage = minCoverage; coverage <= maxCoverage; coverage += coverageStep)
                {
                    // iterate over purity range 
                    for (int percentPurity = minPercentPurity; percentPurity <= maxPercentPurity; percentPurity += 5)
                    {
                        SomaticCaller.CoveragePurityModel model = new SomaticCaller.CoveragePurityModel(somaticCallerParameters.MaximumCopyNumber);
                        model.DiploidCoverage = coverage;
                        model.Purity = percentPurity / 100f;
                        this.ModelDeviation(centroidsMAF, centroidsCoverage, model, usableSegments, bestNumClusters);
                        this.DiploidModelDistance(model, usableSegments, genomeLength);
                        if (model.Deviation < bestDeviation && model.Ploidy < somaticCallerParameters.MaxAllowedPloidy && model.Ploidy > somaticCallerParameters.MinAllowedPloidy)
                        {
                            bestDeviation = model.Deviation;
                        }
                        // exluce models with unrealistic genome ploidies
                        if (model.Ploidy < somaticCallerParameters.MaxAllowedPloidy && model.Ploidy > somaticCallerParameters.MinAllowedPloidy)
                            allModels.Add(model);
                    }
                }
                if (allModels.Count == 0)
                {
                    throw new UncallableDataException(string.Format("Error with CNV detection - unable to find any viable purity/ploidy model."));
                }

                // New logic for model selection:
                // - First, compute the best model deviation.  This establishes a baseline for how large the deviation is allowed to get in 
                //   an acceptable model.  Allow somewhat higher deviation for targeted data, since we see extra noise there.
                // - Review models.  Discard any with unacceptable deviation.  Note the best attainable % copy number 2 and % normal.
                // - For each model, scale PercentNormal to a range of 0..100 where 100 = the best number seen for any acceptable model.  Similarly
                //   for PercentCN2.  And similarly for DeviationScore: BestDeviation=1, WorstAllowedDeviation=0
                // - Choose a model (with acceptable deviation) which maximizes a score of the form:
                //   PercentNormal + a * PercentCN2 + b * DeviationScore
                double worstAllowedDeviation = bestDeviation * somaticCallerParameters.DeviationFactor;
                double bestCN2 = 0;
                double bestCN2Normal = 0;
                double bestDiploidDistance = 0;
                double heterogeneityIndex = 0;

                // derive max values for scaling
                int counter = 0;
                List<double> deviations = new List<double>();
                foreach (SomaticCaller.CoveragePurityModel model in allModels)
                {
                    if (model.Deviation < worstAllowedDeviation) counter++;
                    deviations.Add(model.Deviation);
                }
                deviations.Sort();
                if (counter < somaticCallerParameters.DeviationIndexCutoff)
                {
                    worstAllowedDeviation = deviations[Math.Min(somaticCallerParameters.DeviationIndexCutoff, deviations.Count - 1)];
                }

                double bestAccuracyDeviation = double.MaxValue;
                double bestPrecisionDeviation = double.MaxValue;
                // derive max values for scaling
                foreach (SomaticCaller.CoveragePurityModel model in allModels)
                {
                    bestAccuracyDeviation = Math.Min(bestAccuracyDeviation, model.AccuracyDeviation);
                    bestPrecisionDeviation = Math.Min(bestPrecisionDeviation, model.PrecisionDeviation);
                    if (model.Deviation > worstAllowedDeviation) continue;
                    if (model.PercentCN[2] > bestCN2) bestCN2 = model.PercentCN[2];
                    if (model.PercentNormal > bestCN2Normal) bestCN2Normal = model.PercentNormal;
                    if (model.DiploidDistance > bestDiploidDistance) bestDiploidDistance = model.DiploidDistance;
                }

                // coarse search to find best ploidy and purity model  
                List<SomaticCaller.CoveragePurityModel> bestModels = new List<SomaticCaller.CoveragePurityModel>();
                SomaticCaller.CoveragePurityModel bestModel = null;
                double bestScore = 0;
                // holds scores for all models
                List<double> scores = new List<double>();
                // save all purity and ploidy models to a file 
                string debugPath = Path.Combine(this.TempFolder, "PurityModel.txt");
                using (StreamWriter debugWriter = new StreamWriter(debugPath))
                {
                    debugWriter.Write("#Purity\tDiploidCoverage\t");
                    debugWriter.Write("Deviation\tAccuracyDeviation\tPrecisionDeviation\tWorstAllowedDeviation\tAccDev/best\tPrecDev/best\t");
                    debugWriter.Write("DeviationScore\tScore\tPloidy\t");
                    debugWriter.Write("Normal\tNormal/best\tCN2\tCN2/Best\t");
                    debugWriter.Write("DiploidDistance\tDiploidDistance/Best\t");
                    debugWriter.Write("HeterogeneityIndex\tClusterDeviation");
                    debugWriter.WriteLine();
                    foreach (SomaticCaller.CoveragePurityModel model in allModels)
                    {

                        // Filter models with unacceptable deviation:
                        if (model.Deviation > worstAllowedDeviation) continue;
                        // Transform purity into Weighting Factor to penalize abnormal ploidies at low purity: 
                        // (1.5 - 0.5) = minmax range of the new weighting scale; (1.0 - 0.2) = minmax range of the purity values 
                        // This transformation leads a maximal lowPurityWeightingFactor value of 1.5 for the lowest purity model and a minimal value of 0.75 for the highest purity model 
                        double lowPurityWeightingFactor = 1.5 / ((1.5 - 0.5) / (1.0 - 0.2) * (model.Purity - 0.2) + 1.0);
                        double score = somaticCallerParameters.PercentNormal2WeightingFactor * model.PercentNormal / Math.Max(0.01, bestCN2Normal);
                        if (model.HeterogeneityIndex.HasValue && this.IsEnrichment)
                        {
                            heterogeneityIndex = Math.Max(0.5, Math.Min(4.0, (double)model.HeterogeneityIndex));
                            heterogeneityIndex = Math.Min(1.0, (double)model.HeterogeneityIndex) / Math.Max(1.0, (double)model.HeterogeneityIndex);
                        }

                        score += lowPurityWeightingFactor * somaticCallerParameters.CN2WeightingFactor * model.PercentCN[2] / Math.Max(0.01, bestCN2);
                        score += somaticCallerParameters.DeviationScoreWeightingFactor * (worstAllowedDeviation - model.Deviation) / (worstAllowedDeviation - bestDeviation);
                        score += somaticCallerParameters.DiploidDistanceScoreWeightingFactor * model.DiploidDistance / Math.Max(0.01, bestDiploidDistance);
                        score += somaticCallerParameters.HeterogeneityScoreWeightingFactor * heterogeneityIndex;
                        scores.Add(score);

                        bestModels.Add(model);
                        // write to file
                        debugWriter.Write("{0}\t{1}\t", (int)Math.Round(100 * model.Purity), model.DiploidCoverage);
                        debugWriter.Write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t", model.Deviation, model.AccuracyDeviation, model.PrecisionDeviation,
                            worstAllowedDeviation, model.AccuracyDeviation / bestAccuracyDeviation, model.PrecisionDeviation / bestPrecisionDeviation);
                        debugWriter.Write("{0}\t{1}\t{2}\t", (worstAllowedDeviation - model.Deviation) / (worstAllowedDeviation - bestDeviation),
                            score, model.Ploidy);
                        debugWriter.Write("{0}\t{1}\t{2}\t{3}\t", model.PercentNormal, model.PercentNormal / Math.Max(0.01, bestCN2Normal),
                            model.PercentCN[2], model.PercentCN[2] / Math.Max(0.01, bestCN2));
                        debugWriter.Write("{0}\t{1}\t{2}\t{3}", model.DiploidDistance, model.DiploidDistance / Math.Max(0.01, bestDiploidDistance), heterogeneityIndex, model.ClusterDeviation);
                        debugWriter.WriteLine();

                        if (score > bestScore)
                        {
                            bestModel = model;
                            bestScore = score;
                        }
                    }
                }
                // sort list and return indices
                var sortedScores = scores.Select((x, i) => new KeyValuePair<double, int>(x, i)).OrderBy(x => x.Key).ToList();
                List<double> scoresValue = sortedScores.Select(x => x.Key).ToList();
                List<int> scoresIndex = sortedScores.Select(x => x.Value).ToList();

                // interModelDistance shows genome edit distance between the best model and other top models (defined by MaximumRelatedModels). 
                // The premise is that if the top models provide widely different genome baseline (leading to high interModelDistance), 
                // the overall modeling approach might be more unstable.
                double interModelDistance = 0;
                // start at one since model #0 is the highest scoring model to compare to
                for (int i = 1; i < Math.Min(scoresIndex.Count, SomaticCallerParameters.MaximumRelatedModels); i++)
                {
                    interModelDistance += CalculateModelDistance(bestModels[scoresIndex[0]], bestModels[scoresIndex[i]], usableSegments, genomeLength);
                }
                interModelDistance /= (double)SomaticCallerParameters.MaximumRelatedModels;

                Console.WriteLine(">>> Initial model: Deviation {0:F5}, coverage {1}, purity {2:F1}%, CN2 {3:F2}", bestModel.Deviation,
                        bestModel.DiploidCoverage, 100 * bestModel.Purity, bestModel.PercentCN[2]);

                // Refine search: Smaller step sizes in the neighborhood of the initial model.
                if (this.userPloidy != null)
                {
                    minCoverage = maxCoverage = (int)GetDiploidCoverage(medianCoverageLevel, this.userPloidy.Value);
                }
                else
                {
                    minCoverage = (int)Math.Round(bestModel.DiploidCoverage) - 5;
                    maxCoverage = (int)Math.Round(bestModel.DiploidCoverage) + 5;
                }
                if (this.userPurity != null)
                {
                    minPercentPurity = maxPercentPurity = (int)(this.userPurity.Value * 100);
                }
                else
                {
                    minPercentPurity = Math.Max(20, (int)Math.Round(bestModel.Purity * 100) - 10);
                    maxPercentPurity = Math.Min(100, (int)Math.Round(bestModel.Purity * 100) + 10); // %%% magic numbers
                }
                bestDeviation = double.MaxValue;
                bestModel = null;
                for (int coverage = minCoverage; coverage <= maxCoverage; coverage++)
                {
                    for (int percentPurity = minPercentPurity; percentPurity <= maxPercentPurity; percentPurity++)
                    {
                        SomaticCaller.CoveragePurityModel model = new SomaticCaller.CoveragePurityModel(somaticCallerParameters.MaximumCopyNumber);
                        model.DiploidCoverage = coverage;
                        model.Purity = percentPurity / 100f;
                        this.ModelDeviation(centroidsMAF, centroidsCoverage, model, usableSegments, bestNumClusters);
                        if (bestModel == null || model.Deviation < bestModel.Deviation)
                        {
                            bestModel = model;
                        }
                    }
                }
                string debugPathClusterModel = Path.Combine(this.TempFolder, "ClusteringModel.txt");
                string debugPathCNVModeling = Path.Combine(this.TempFolder, "CNVModeling.txt");

                ModelDeviation(centroidsMAF, centroidsCoverage, bestModel, usableSegments, bestNumClusters, debugPathClusterModel, true, debugPathCNVModeling);
                Console.WriteLine();
                Console.WriteLine(">>> Refined model: Deviation {0:F5}, coverage {1}, purity {2:F1}%", bestModel.Deviation,
                    bestModel.DiploidCoverage, bestModel.Purity * 100);
                Console.WriteLine();


                {
                    foreach (SegmentPloidy ploidy in AllPloidies)
                    {
                        ploidy.Omega = 0;
                        ploidy.Mu = null;
                        ploidy.Sigma = null;
                    }
                }
                if (!bestModel.InterModelDistance.HasValue)
                {
                    bestModel.InterModelDistance = interModelDistance;
                }
                return bestModel;
            }
        }

        private static double GetDiploidCoverage(int medianCoverageLevel, float ploidy)
        {
            return medianCoverageLevel / ploidy * 2.0;
        }

        /// <summary>
        /// Assign a SegmentPloidy to each CanvasSegment, based on which model matches this segment best:
        /// </summary>
        protected void AssignPloidyCalls()
        {

            foreach (CanvasSegment segment in this.Segments)
            {
                // Compute (MAF, Coverage) for this segment:
                List<double> MAF = new List<double>();
                foreach (float VF in segment.VariantFrequencies) MAF.Add(VF > 0.5 ? 1 - VF : VF);
                double medianCoverage = CanvasCommon.Utilities.Median(segment.Counts);
                MAF.Sort();

                // For segments with (almost) no variants alleles at all, we'll assign them a dummy MAF, and 
                // we simply won't consider MAF when determining the closest ploidy:
                double? medianMAF = null;

                if (MAF.Count >= 10)
                {
                    medianMAF = MAF[MAF.Count / 2];
                }

                // Find the closest ploidy, using euclidean distance:
                double bestDistance = double.MaxValue;
                double secondBestDistance = double.MaxValue;
                SegmentPloidy bestPloidy = null;
                SegmentPloidy secondBestPloidy = null;

                foreach (SegmentPloidy ploidy in AllPloidies)
                {
                    double distance = GetModelDistance(ploidy.MixedCoverage, medianCoverage, medianMAF, ploidy.MixedMinorAlleleFrequency);

                    if (distance < bestDistance)
                    {
                        secondBestDistance = bestDistance;
                        secondBestPloidy = bestPloidy;
                        bestDistance = distance;
                        bestPloidy = ploidy;
                    }
                    else if (distance < secondBestDistance)
                    {
                        secondBestDistance = distance;
                        secondBestPloidy = ploidy;
                    }
                }
                segment.CopyNumber = bestPloidy.CopyNumber;
                segment.SecondBestCopyNumber = secondBestPloidy.CopyNumber;

                // If we don't have variant frequencies, then don't report any major chromosome count
                segment.MajorChromosomeCount = null;
                if (medianMAF.HasValue)
                {
                    segment.MajorChromosomeCount = bestPloidy.MajorChromosomeCount;
                }
                segment.ModelDistance = bestDistance;
                segment.RunnerUpModelDistance = secondBestDistance;

                // Special logic for extra-high coverage: Adjust our copy number call to be scaled by coverage,
                // and adjust our model distance as well.
                if (segment.CopyNumber == somaticCallerParameters.MaximumCopyNumber)
                {
                    double coverageRatio = segment.MeanCount / Model.DiploidCoverage;
                    int referenceCN = 2;
                    if (this.ReferencePloidy != null) referenceCN = this.ReferencePloidy.GetReferenceCopyNumber(segment);
                    // Formula for reference:
                    // 2*coverageRatio = ObservedCopyNumber = Purity*TumorCopyNumber + (1-Purity)*ReferenceCopyNumber
                    double estimateCopyNumber = (2 * coverageRatio - referenceCN * (1 - Model.Purity)) / Model.Purity;
                    int estimateCN = (int)Math.Round(estimateCopyNumber);
                    if (estimateCN > somaticCallerParameters.MaximumCopyNumber)
                    {
                        segment.CopyNumber = estimateCN;
                        segment.MajorChromosomeCount = null;
                        double coverage = Model.DiploidCoverage * ((1 - Model.Purity) + Model.Purity * estimateCN / 2f);
                        segment.ModelDistance = Math.Abs(segment.MeanCount - coverage) * this.CoverageWeightingFactor;
                    }
                }

            }
        }

        /// <summary>
        /// Adjust CN in segments from heterogeneous clusters when the following conditions are satisfied:
        ///      1.Best CN = 2, but second best ploidy suggest CN = 1 or 3
        ///      2.Deviations of best and second best ploidy model are very similar
        ///      3.Purity is reasonably high
        ///      4.There’s an evidence that tumor is heterogeneous
        /// In such a case CN = 2 most likely represents are heterogeneous variant
        /// </summary>
        protected void AdjustPloidyCalls()
        {


            Console.WriteLine();
            Console.WriteLine(">>> HeterogeneousSegmentsSignature: {0}", HeterogeneousSegmentsSignature.Count);
            Console.WriteLine();

            foreach (CanvasSegment segment in this.Segments)
            {
                if (this.HeterogeneousSegmentsSignature.BinarySearch(segment.Begin + segment.End + segment.Counts.Count) >= 0 && Model.Purity > 0.2f && Math.Max(segment.ModelDistance, Double.MinValue) / segment.RunnerUpModelDistance > somaticCallerParameters.DistanceRatio)
                {
                    if (segment.CopyNumber == 2 && (segment.SecondBestCopyNumber == 1 || segment.SecondBestCopyNumber == 3))
                    {

                        int tmpCopyNumber = segment.SecondBestCopyNumber;
                        segment.SecondBestCopyNumber = segment.CopyNumber; // indicator that CNs have swapped
                        segment.CopyNumber = tmpCopyNumber;
                        segment.cnSwaped = "Y";

                    }
                }
            }
        }

        /// <summary>
        /// Assign a SegmentPloidy to each CanvasSegment, based on which model matches this segment best:
        /// </summary>
        protected void AssignPloidyCallsGaussianMixture()
        {
            // For segments with (almost) no variants alleles at all, we'll assign them a dummy MAF, and 
            // we simply won't consider MAF when determining the closest ploidy:
            double dummyMAF = -1;

            foreach (CanvasSegment segment in this.Segments)
            {
                // Compute (MAF, Coverage) for this segment:
                List<double> MAF = new List<double>();
                foreach (float VF in segment.VariantFrequencies) MAF.Add(VF > 0.5 ? 1 - VF : VF);
                double medianCoverage = CanvasCommon.Utilities.Median(segment.Counts);
                MAF.Sort();
                double medianMAF = dummyMAF;
                if (MAF.Count >= 10)
                {
                    medianMAF = MAF[MAF.Count / 2];
                }

                Dictionary<SegmentPloidy, double> posteriorProbabilities = GaussianMixtureModel.EMComputePosteriorProbs(AllPloidies, medianMAF, medianCoverage);

                // Find the closest ploidy. 
                double bestProbability = 0;
                SegmentPloidy bestPloidy = null;
                foreach (SegmentPloidy ploidy in AllPloidies)
                {
                    if (bestPloidy == null || posteriorProbabilities[ploidy] > bestProbability)
                    {
                        bestProbability = posteriorProbabilities[ploidy];
                        bestPloidy = ploidy;
                    }
                }

                // Sanity-check: If we didn't find anything with probability > 1, then fall back to the simplest possible
                // thing: Call purely on coverage.
                if (bestProbability == 0)
                {
                    segment.CopyNumber = (int)Math.Round(2 * medianCoverage / Model.DiploidCoverage);
                    segment.MajorChromosomeCount = segment.CopyNumber / 2;
                }
                else
                {
                    segment.CopyNumber = bestPloidy.CopyNumber;
                    segment.MajorChromosomeCount = bestPloidy.MajorChromosomeCount;
                }
            }
        }

        /// <summary>
        /// Assign copy number calls to segments.  And, produce extra headers for the CNV vcf file, giving the 
        /// overall estimated purity and ploidy.
        /// </summary>
        protected List<string> CallCNVUsingSNVFrequency(double? localSDmertic, string referenceFolder, CanvasSomaticClusteringMode clusteringMode)
        {
            List<string> Headers = new List<string>();
            if (this.CNOracle != null)
            {
                this.DerivePurityEstimateFromVF();
            }

            // Get genome length.
            GenomeMetadata genomeMetaData = null;
            genomeMetaData = new GenomeMetadata();
            genomeMetaData.Deserialize(Path.Combine(referenceFolder, "GenomeSize.xml"));

            // Derive a model of diploid coverage, and overall tumor purity:
            this.Model = ModelOverallCoverageAndPurity(genomeMetaData.Length, clusteringMode);

            // Make preliminary ploidy calls for all segments.  For those segments which fit their ploidy reasonably well,
            // accumulate information about the MAF by site and coverage by bin.  
            this.HeterogeneousSegmentsSignature.Sort();

            if (AllPloidies.First().Sigma == null)
            {
                AssignPloidyCalls();
                // Do not run heterogeneity adjustment on enrichment data
                if (!this.IsEnrichment)
                    AdjustPloidyCalls();
            }
            else
            {
                AssignPloidyCallsGaussianMixture();
            }

            // If the somatic SNV/indel file was provided, then we use it to derive another estimate of purity.
            // And, if we didn't make many CNV calls, then we report this estimate, instead of the estimate derived from
            // our overall model.
            if (!string.IsNullOrEmpty(SomaticVCFPath))
            {
                try
                {
                    double SNVPurityEstimate = EstimatePurityFromSomaticSNVs();
                    this.SelectPurityEstimate(SNVPurityEstimate, genomeMetaData.Length);
                }
                catch (Exception e)
                {
                    Console.Error.WriteLine("* Error deriving purity estimate from somatic SNVs.  Details:\n{0}", e.ToString());
                }
            }

            // Add some extra information to the vcf file header:
            Headers.Add(string.Format("##EstimatedTumorPurity={0:F2}", this.Model.Purity));
            Headers.Add(string.Format("##PurityModelFit={0:F4}", this.Model.Deviation));
            Headers.Add(string.Format("##InterModelDistance={0:F4}", this.Model.InterModelDistance));
            Headers.Add(string.Format("##LocalSDmetric={0:F2}", localSDmertic));
            Headers.Add(string.Format("##Heterogeneity={0:F2}", this.Model.HeterogeneityIndex));
            return Headers;
        }

        /// <summary>
        /// Provide an estimated chromosome count: Sum up the total number of autosomes and allosomes after taking into 
        /// account all the PF CNV calls.
        /// </summary>
        private double EstimateChromosomeCount()
        {
            int[] baseCountByCopyNumber = new int[somaticCallerParameters.MaximumCopyNumber + 1];
            string currentChromosome = null;
            double overallCount = 0;
            foreach (CanvasSegment segment in this.Segments)
            {
                if (segment.Chr != currentChromosome)
                {
                    if (currentChromosome != null)
                        overallCount += GetWeightedChromosomeCount(baseCountByCopyNumber);
                    Array.Clear(baseCountByCopyNumber, 0, baseCountByCopyNumber.Length);
                    currentChromosome = segment.Chr;
                }
                if (segment.Filter != "PASS") continue;
                baseCountByCopyNumber[Math.Min(segment.CopyNumber, somaticCallerParameters.MaximumCopyNumber)] += (segment.End - segment.Begin);
            }
            overallCount += GetWeightedChromosomeCount(baseCountByCopyNumber);
            return overallCount;
        }

        private double GetWeightedChromosomeCount(int[] baseCountByCopyNumber)
        {
            double weightedMean = 0;
            double weight = 0;
            for (int CN = 0; CN < baseCountByCopyNumber.Length; CN++)
            {
                weight += baseCountByCopyNumber[CN];
                weightedMean += baseCountByCopyNumber[CN] * CN;
            }
            if (weight == 0) return 0;
            return weightedMean / weight;
        }


        /// <summary>
        /// Override our model's purity estimate with the estimate from somatic SNVs, if we didn't call enough CNVs
        /// to get a reasonable estimate of purity from the CNV data:
        /// </summary>
        protected void SelectPurityEstimate(double SNVPurityEstimate, long genomeLength)
        {
            Console.WriteLine("Purity estimates: {0:F4} from CNVs, {1:F4} from SNVs", this.Model.Purity, SNVPurityEstimate);
            double fractionAbnormal = 0;
            double totalWeight = 0;
            foreach (CanvasSegment segment in this.Segments)
            {
                totalWeight += segment.End - segment.Begin;
                if (segment.CopyNumber != 2 || segment.MajorChromosomeCount != 1)
                {
                    fractionAbnormal += segment.End - segment.Begin;
                }
            }
            fractionAbnormal /= totalWeight;
            Console.WriteLine(">>>Fraction abnormal CNV is {0:F4}", fractionAbnormal);
            if (fractionAbnormal < 0.07 && !double.IsNaN(SNVPurityEstimate) && this.Model.Purity < 0.5)
            {
                Console.WriteLine(">>>Override purity estimate to {0:F4}", SNVPurityEstimate);
                this.Model.Purity = SNVPurityEstimate;
            }
        }

        /// <summary>
        /// Load the somatic SNV output (from Strelka), and use it to dervie an estimate of overall purity.
        /// Reference: https://ukch-confluence.illumina.com/display/collab/Estimate+purity+of+tumour+samples+using+somatic+SNVs
        /// </summary>
        protected double EstimatePurityFromSomaticSNVs()
        {
            Dictionary<string, List<CanvasSegment>> segmentsByChromosome = CanvasSegment.GetSegmentsByChromosome(Segments);
            int recordCount = 0;
            List<float> variantFrequencies = new List<float>();
            using (VcfReader reader = new VcfReader(this.SomaticVCFPath))
            {
                foreach (VcfVariant variant in reader.GetVariants())
                {
                    recordCount++;
                    // Skip non-PF variants:
                    if (variant.Filters != "PASS") continue;
                    // Skip everything but SNVs:
                    if (variant.ReferenceAllele.Length > 1 || variant.VariantAlleles.Length != 1 || variant.VariantAlleles[0].Length != 1
                        || variant.VariantAlleles[0] == ".")
                    {
                        continue;
                    }

                    string refTagName = string.Format("{0}U", variant.ReferenceAllele[0]);
                    string altTagName = string.Format("{0}U", variant.VariantAlleles[0][0]);
                    string[] counts = variant.GenotypeColumns.Last()[refTagName].Split(',');
                    int refCount = 0;
                    foreach (string bit in counts) refCount += int.Parse(bit);
                    int altCount = 0;
                    counts = variant.GenotypeColumns.Last()[altTagName].Split(',');
                    foreach (string bit in counts) altCount += int.Parse(bit);
                    float VF = altCount / (float)(altCount + refCount);
                    if (VF >= 0.5) continue;
                    variantFrequencies.Add(VF);
                }
            }
            Console.WriteLine(">>>Loaded {0} somatic variants; saved {1} somatic SNV frequencies", recordCount, variantFrequencies.Count);
            // Don't bother providing an estimate if we have very few events:
            if (variantFrequencies.Count < 100)
            {
                return double.NaN;
            }
            double mean = variantFrequencies.Average();
            double estimatedPurity = Math.Min(1, mean * 2);
            Console.WriteLine(">>>Estimated tumor purity of {0} from somatic SNV calls", estimatedPurity);

            return estimatedPurity;
        }

        /// <summary>
        /// Generate a table listing segments (and several features), and noting which are accurate (copy number 
        /// exactly matches truth set) or directionally accurate (copy number and truth set are both <2, both =2, or both >2)
        /// This table will become our collection of feature vectors for training q-scores!
        /// </summary>
        private void GenerateReportVersusKnownCN()
        {
            string debugPath = Path.Combine(this.TempFolder, "CallsVersusKnownCN.txt");
            using (StreamWriter writer = new StreamWriter(debugPath))
            {
                writer.Write("#Chr\tBegin\tEnd\tTruthSetCN\tTruthSetClonality\tBestCN\tSecondBestCN\tisCNswapped\t");
                writer.Write("Accurate\tDirectionAccurate\t");
                writer.Write("LogLength\tLogBinCount\tBinCount\tBinCV\tModelDistance\tRunnerUpModelDistance\t");
                writer.Write("MafCount\tMafMean\tMafCv\tLogMafCv\tCopyNumber\tMCC\t");
                writer.Write("DistanceRatio\tLogMafCount\t");
                writer.Write("IsHeterogeneous\tModelPurity\tModelDeviation\t");
                writer.Write("QScoreLinearFit\tQScoreGeneralizedLinearFit\tQScoreLogistic\t");
                writer.WriteLine();
                foreach (CanvasSegment segment in this.Segments)
                {
                    int CN = this.GetKnownCNForSegment(segment);
                    double Heterogeneity = this.GetKnownClonalityForSegment(segment);
                    if (CN < 0) continue;
                    if (segment.End - segment.Begin < 5000) continue;
                    List<float> MAF = new List<float>();
                    foreach (float VF in segment.VariantFrequencies)
                    {
                        MAF.Add(VF > 0.5 ? 1 - VF : VF);
                    }
                    MAF.Sort();
                    float MedianMAF = -1;
                    if (MAF.Count > 0)
                        MedianMAF = MAF[MAF.Count / 2];
                    double medianCoverage = CanvasCommon.Utilities.Median(segment.Counts);
                    string accurateFlag = "N";
                    if (CN == segment.CopyNumber) accurateFlag = "Y";
                    string directionAccurateFlag = "N";
                    if ((CN < 2 && segment.CopyNumber < 2) ||
                        (CN == 2 && segment.CopyNumber == 2) ||
                        (CN > 2 && segment.CopyNumber > 2))
                        directionAccurateFlag = "Y";
                    writer.Write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t", segment.Chr, segment.Begin, segment.End, CN, Heterogeneity, segment.CopyNumber, segment.SecondBestCopyNumber, segment.cnSwaped);
                    writer.Write("{0}\t{1}\t", accurateFlag, directionAccurateFlag);
                    writer.Write("{0}\t", Math.Log(segment.End - segment.Begin));
                    writer.Write("{0}\t", segment.GetQScorePredictor(CanvasSegment.QScorePredictor.LogBinCount));
                    writer.Write("{0}\t", segment.GetQScorePredictor(CanvasSegment.QScorePredictor.BinCount));
                    writer.Write("{0}\t", segment.GetQScorePredictor(CanvasSegment.QScorePredictor.BinCv));
                    writer.Write("{0}\t", segment.GetQScorePredictor(CanvasSegment.QScorePredictor.ModelDistance));
                    writer.Write("{0}\t", segment.GetQScorePredictor(CanvasSegment.QScorePredictor.RunnerUpModelDistance));
                    writer.Write("{0}\t", segment.GetQScorePredictor(CanvasSegment.QScorePredictor.MafCount));
                    writer.Write("{0}\t", segment.GetQScorePredictor(CanvasSegment.QScorePredictor.MafMean));
                    writer.Write("{0}\t", segment.GetQScorePredictor(CanvasSegment.QScorePredictor.MafCv));
                    writer.Write("{0}\t", segment.GetQScorePredictor(CanvasSegment.QScorePredictor.LogMafCv));
                    writer.Write("{0}\t", segment.GetQScorePredictor(CanvasSegment.QScorePredictor.CopyNumber));
                    writer.Write("{0}\t", segment.GetQScorePredictor(CanvasSegment.QScorePredictor.MajorChromosomeCount));
                    writer.Write("{0}\t", segment.GetQScorePredictor(CanvasSegment.QScorePredictor.DistanceRatio));
                    writer.Write("{0}\t", Math.Log(segment.GetQScorePredictor(CanvasSegment.QScorePredictor.MafCount)));
                    if (this.HeterogeneousSegmentsSignature.BinarySearch(segment.Begin + segment.End + segment.Counts.Count) >= 0)
                        writer.Write("{0}\t", "Y");
                    else
                    {
                        writer.Write("{0}\t", "N");
                    }
                    writer.Write("{0}\t", Model.Purity);
                    writer.Write("{0}\t", Model.Deviation);
                    double score = segment.ComputeQScore(CanvasSegment.QScoreMethod.BinCountLinearFit);
                    writer.Write("{0}\t", score);
                    score = segment.ComputeQScore(CanvasSegment.QScoreMethod.GeneralizedLinearFit);
                    writer.Write("{0}\t", score);
                    score = segment.ComputeQScore(CanvasSegment.QScoreMethod.Logistic);
                    writer.Write("{0}\t", score);
                    writer.WriteLine();
                }
            }
            Console.WriteLine(">>> Wrote report of CNV calls versus reference calls to {0}", debugPath);
        }

        /// <summary>
        /// Developer debug method:
        /// - Split each truth interval to have at least the same segmentation as the called segments
        ///   (Note: We are intentionally ignoring segments - or parts thereof - called in areas not defined in the Truth set)
        /// - For each QScore method:
        ///   - Report these new intervals and associated QScores to an extended report output file
        ///   - Generate ROC output
        /// </summary>
        private void GenerateExtendedReportVersusKnownCN()
        {
            Dictionary<string, List<CNInterval>> resegmentedKnownCN = new Dictionary<string, List<CNInterval>>();

            // Copy KnownCN entries to working container
            foreach (string chr in this.CNOracle.KnownCN.Keys)
            {
                resegmentedKnownCN[chr] = new List<CNInterval>();
                foreach (CNInterval interval in this.CNOracle.KnownCN[chr])
                {
                    CNInterval newInterval = new CNInterval();
                    newInterval.Start = interval.Start;
                    newInterval.End = interval.End;
                    newInterval.CN = interval.CN;
                    resegmentedKnownCN[chr].Add(newInterval);
                }
            }

            // Split each truth interval to match the segments' breakpoints
            foreach (CanvasSegment segment in this.Segments)
            {
                if (!resegmentedKnownCN.ContainsKey(segment.Chr)) continue;
                for (int i = 0; i < resegmentedKnownCN[segment.Chr].Count; i++) // Using for loop instead of foreach because we add items to the list
                {
                    CNInterval interval = resegmentedKnownCN[segment.Chr][i];
                    if (interval.Start == segment.Begin && interval.End == segment.End) break; // perfect segment-knownCN match
                    if (interval.Start >= segment.End || interval.End <= segment.Begin) continue; // segment completely outside this knownCN

                    // If necessary, split interval at segment.Begin position (i.e. extract sub-interval preceding segment)
                    if (segment.Begin > interval.Start)
                    {
                        CNInterval newInterval = new CNInterval();
                        newInterval.Start = interval.Start;
                        newInterval.End = segment.Begin;
                        newInterval.CN = interval.CN;
                        interval.Start = newInterval.End;
                        resegmentedKnownCN[segment.Chr].Add(newInterval);
                    }

                    // If necessary, split interval at segment.End position (i.e. extract sub-interval following segment)
                    if (segment.End < interval.End)
                    {
                        CNInterval newInterval = new CNInterval();
                        newInterval.Start = segment.End;
                        newInterval.End = interval.End;
                        newInterval.CN = interval.CN;
                        interval.End = newInterval.Start;
                        resegmentedKnownCN[segment.Chr].Add(newInterval);
                    }
                }
            }

            // Sort list of new intervals by starting position, just for prettiness
            foreach (List<CNInterval> list in resegmentedKnownCN.Values)
            {
                list.Sort((i1, i2) => i1.Start.CompareTo(i2.Start));
            }

            // Generate ROC output data for each QScore method
            foreach (CanvasSegment.QScoreMethod qscoreMethod in CanvasSegment.QScoreMethod.GetValues(typeof(CanvasSegment.QScoreMethod)))
            {
                GenerateReportAndRocDataForQscoreMethod(qscoreMethod, resegmentedKnownCN);
            }
        }

        /// <summary>
        /// Developer debug method: ROC curve data generation
        /// - Report all intervals, associated QScores and QScore predictor values to an extended report output file
        /// - Report called (i.e. TP+FP) intervals grouped by QScore
        /// - Generate 2 ROC outputs
        ///   - ROC_intervals: FP vs TP rate, unit=1 interval
        ///   - ROC_bases:     FP vs TP rate, unit=1 base
        ///   (Note: In both cases, we ignore intervals shorter than 1kb as most of them are due to imprecise ends of segments, which we don't want to give any weight to)
        /// </summary>
        private void GenerateReportAndRocDataForQscoreMethod(CanvasSegment.QScoreMethod qscoreMethod, Dictionary<string, List<CNInterval>> resegmentedKnownCN)
        {
            // Create map interval->{segment+qscore}, ignoring intervals shorter than 1kb
            Dictionary<CNInterval, Tuple<CanvasSegment, int>> Interval2Segment = new Dictionary<CNInterval, Tuple<CanvasSegment, int>>();
            foreach (string chr in resegmentedKnownCN.Keys)
            {
                foreach (CNInterval interval in resegmentedKnownCN[chr])
                {
                    foreach (CanvasSegment segment in this.Segments)
                    {
                        if (segment.Chr == chr && (segment.Begin == interval.Start || segment.End == interval.End))
                        {
                            if (interval.End - interval.Start >= 1000)
                                Interval2Segment[interval] = new Tuple<CanvasSegment, int>(segment, segment.ComputeQScore(qscoreMethod));
                        }
                    }
                }
            }

            // Classify intervals by QScore
            List<List<CNInterval>> intervalsByQScore = new List<List<CNInterval>>();
            foreach (CNInterval interval in Interval2Segment.Keys)
            {
                int qscore = Interval2Segment[interval].Item2;
                // Resize list to hold this qscore's entries
                while (qscore >= intervalsByQScore.Count())
                {
                    intervalsByQScore.Add(new List<CNInterval>());
                }
                intervalsByQScore[qscore].Add(interval);
            }

            // Output data as ExtendedCallsVersusKnownCN.txt
            string debugPath = Path.Combine(this.OutputFolder, "qscore_" + qscoreMethod.ToString() + "_ExtendedCallsVersusKnownCN.txt");
            using (StreamWriter writer = new StreamWriter(debugPath))
            {
                writer.Write("#Chr\tBegin\tEnd\tTruthSetCN\tCalledCN\tMajorChromCount\tQScore\tInfo");
                foreach (CanvasSegment.QScorePredictor predictorId in CanvasSegment.QScorePredictor.GetValues(typeof(CanvasSegment.QScorePredictor)))
                {
                    writer.Write("\tPredictor_{0}", predictorId.ToString());
                }
                writer.WriteLine("");

                foreach (string chr in resegmentedKnownCN.Keys)
                {
                    foreach (CNInterval interval in resegmentedKnownCN[chr])
                    {
                        if (Interval2Segment.ContainsKey(interval))
                        {
                            CanvasSegment segment = Interval2Segment[interval].Item1;
                            int qscore = Interval2Segment[interval].Item2;
                            string info = (interval.CN == segment.CopyNumber) ? "Correct" : "Incorrect";
                            writer.Write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}", chr, interval.Start, interval.End, interval.CN, segment.CopyNumber, segment.MajorChromosomeCount, qscore, info);
                            foreach (CanvasSegment.QScorePredictor predictorId in CanvasSegment.QScorePredictor.GetValues(typeof(CanvasSegment.QScorePredictor)))
                            {
                                writer.Write("\t{0}", segment.GetQScorePredictor(predictorId));
                            }
                            writer.WriteLine("");
                        }
                        else
                        {
                            string info = "Missing";
                            int CN = -1;
                            int majorChromosomeCount = -1;
                            writer.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", chr, interval.Start, interval.End, interval.CN, CN, majorChromosomeCount, info);
                        }
                    }
                }
            }

            // Output data by QScore
            debugPath = Path.Combine(this.OutputFolder, "qscore_" + qscoreMethod + "_cnaPerQscore.txt");
            using (StreamWriter writer = new StreamWriter(debugPath))
            {
                writer.WriteLine("#Chr\tBegin\tEnd\tTruthSetCN\tCalledCN\tMajorChromCount\tMedianMAF\tMedianCoverage\tQScore\tInfo");
                for (int qscore = 0; qscore < intervalsByQScore.Count(); qscore++)
                {
                    foreach (CNInterval interval in intervalsByQScore[qscore])
                    {
                        CanvasSegment segment = Interval2Segment[interval].Item1;
                        string info = (interval.CN == segment.CopyNumber) ? "Correct" : "Incorrect";
                        writer.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}", segment.Chr, interval.Start, interval.End, interval.CN, segment.CopyNumber, segment.MajorChromosomeCount, qscore, info);
                    }
                }
            }

            // ROC output per interval
            debugPath = Path.Combine(this.OutputFolder, "qscore_" + qscoreMethod + "_ROC_intervals.txt");
            GenerateRocOutput(debugPath, intervalsByQScore, Interval2Segment, false, false);

            // ROC output per base
            debugPath = Path.Combine(this.OutputFolder, "qscore_" + qscoreMethod + "_ROC_bases.txt");
            GenerateRocOutput(debugPath, intervalsByQScore, Interval2Segment, true, false);
        }

        /// <summary>
        /// Developer debug method: ROC data output
        /// </summary>
        private void GenerateRocOutput(string debugPath, List<List<CNInterval>> intervalsByQScore, Dictionary<CNInterval, Tuple<CanvasSegment, int>> Interval2Segment, bool countBases, bool scaleUsingMaxFp)
        {
            // Gather some overall counts
            long totalFP = 0;
            long totalTP = 0;
            foreach (CNInterval interval in Interval2Segment.Keys)
            {
                CanvasSegment segment = Interval2Segment[interval].Item1;
                if (interval.CN == segment.CopyNumber)
                    totalTP += countBases ? interval.End - interval.Start : 1;
                else
                    totalFP += countBases ? interval.End - interval.Start : 1;
            }

            using (StreamWriter writer = new StreamWriter(debugPath))
            {
                writer.WriteLine("#QScore\tFP\tTP\tFPR(retained FP)\tTPR");
                long TP = 0;
                long FP = 0;
                long maxFP = 0;
                double lastFPR = 0;
                double lastTPR = 0;
                double rocArea = 0;
                for (int qscore = intervalsByQScore.Count() - 1; qscore >= 0; qscore--)
                {
                    foreach (CNInterval interval in intervalsByQScore[qscore])
                    {
                        CanvasSegment segment = Interval2Segment[interval].Item1;
                        if (interval.CN == segment.CopyNumber)
                            TP += countBases ? interval.End - interval.Start : 1;
                        else
                            FP += countBases ? interval.End - interval.Start : 1;
                    }
                    double FPR;
                    if (scaleUsingMaxFp)
                    {
                        if (qscore == intervalsByQScore.Count() - 1)
                            maxFP = FP;
                        FPR = ((double)FP - maxFP) / (totalFP - maxFP);
                    }
                    else
                        FPR = ((double)FP) / totalFP;
                    double TPR = ((double)TP) / totalTP;
                    writer.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", qscore, FP, TP, FPR, TPR);

                    // ROC Area calculation
                    rocArea += (TPR + lastTPR) / 2 * (FPR - lastFPR);
                    lastFPR = FPR;
                    lastTPR = TPR;
                }
                writer.WriteLine("#Area=\t{0}", rocArea);
            }
        }

        /// <summary>
        /// Developer debug method:
        /// - Loop over intervals where copy number is known
        /// - Intersect these with our segments
        /// - Count the # of bases where we have a correct or incorrect CN call.
        /// - Give "partial credit" if we've got the correct state of <2, =2, >2
        /// </summary>
        private void DebugEvaluateCopyNumberCallAccuracy()
        {
            long correctBases = 0;
            long incorrectBases = 0;
            long allSegmentBases = 0;
            long allTruthBases = 0;
            long partialCreditBases = 0;
            foreach (List<CNInterval> list in this.CNOracle.KnownCN.Values)
            {
                foreach (CNInterval interval in list)
                {
                    allTruthBases += interval.End - interval.Start;
                }
            }
            long[,] ConfusionMatrix = new long[8, 8];
            foreach (CanvasSegment segment in this.Segments)
            {
                allSegmentBases += segment.End - segment.Begin;
                if (!this.CNOracle.KnownCN.ContainsKey(segment.Chr)) continue;
                foreach (CNInterval interval in this.CNOracle.KnownCN[segment.Chr])
                {
                    if (interval.Start > segment.End) continue;
                    if (interval.End < segment.Begin) continue;
                    int start = Math.Max(interval.Start, segment.Begin);
                    int end = Math.Min(interval.End, segment.End);
                    if (interval.CN == segment.CopyNumber)
                    {
                        correctBases += (end - start);
                        if (interval.CN < 8)
                            ConfusionMatrix[interval.CN, segment.CopyNumber] += (end - start);
                    }
                    else
                    {
                        if (interval.CN < 8 && segment.CopyNumber < 8)
                            ConfusionMatrix[interval.CN, segment.CopyNumber] += (end - start);
                        incorrectBases += (end - start);
                    }
                    if ((interval.CN == 2 && segment.CopyNumber == 2) ||
                        (interval.CN < 2 && segment.CopyNumber < 2) ||
                        (interval.CN > 2 && segment.CopyNumber > 2))
                        partialCreditBases += (end - start);
                }
            }
            long totalIntersect = correctBases + incorrectBases;
            Console.WriteLine("Correct bases: {0} Incorrect bases: {1} Accuracy: {2:F2}%, Partial accuracy {3:F2}%", correctBases,
                incorrectBases, 100 * correctBases / (double)(totalIntersect),
                100 * partialCreditBases / (double)totalIntersect);
            Console.WriteLine("Truth set and call set overlap: {0} bases", totalIntersect);
            Console.WriteLine("Intersect {0:F2}% of {1} called bases", 100 * totalIntersect / (double)allSegmentBases, allSegmentBases);
            Console.WriteLine("Intersect {0:F2}% of {1} truth-set bases bases", 100 * totalIntersect / (double)allTruthBases, allTruthBases);
            string debugPath = Path.Combine(this.TempFolder, "CNVConfusionMatrix.txt");
            using (StreamWriter writer = new StreamWriter(debugPath))
            {
                writer.WriteLine("#OverallAccuracy\t{0:F2}\t{1:F2}", 100 * correctBases / (double)totalIntersect,
                    100 * partialCreditBases / (double)totalIntersect);
                // Overall confusion matrix:
                writer.WriteLine("Overall confusion matrix:");
                writer.WriteLine("#Actual \\ Called\t0\t1\t2\t3\t4\t5\t6\t7\t");
                for (int index = 0; index < 8; index++)
                {
                    writer.Write("{0}\t", index);
                    for (int index2 = 0; index2 < 8; index2++)
                        writer.Write("{0}\t", ConfusionMatrix[index, index2]);
                    writer.WriteLine();
                }
                writer.WriteLine();

                writer.WriteLine("Confusion matrix, normalized by ACTUAL copy number:");
                writer.WriteLine("#Actual \\ Called\t0\t1\t2\t3\t4\t5\t6\t7\t");
                for (int index = 0; index < 8; index++)
                {
                    writer.Write("{0}\t", index);
                    float count = 0;
                    for (int index2 = 0; index2 < 8; index2++)
                        count += ConfusionMatrix[index, index2];
                    if (count == 0) count = 1;
                    for (int index2 = 0; index2 < 8; index2++)
                        writer.Write("{0}\t", 100 * ConfusionMatrix[index, index2] / count);
                    writer.WriteLine();
                }
            }
            Console.WriteLine("Confusion matrix written to {0}.", debugPath);
        }

        /// <summary>
        /// Check whether we know the CN for this segment.  Look for a known-CN interval that 
        /// covers (at least half of) this segment.  Return -1 if we don't know its CN.
        /// </summary>
        protected int GetKnownCNForSegment(CanvasSegment segment)
        {
            if (this.CNOracle == null) return -1;
            return this.CNOracle.GetKnownCNForSegment(segment);
        }

        protected double GetKnownClonalityForSegment(CanvasSegment segment)
        {
            if (this.CNOracle == null) return -1;
            return this.CNOracle.GetKnownClonalityForSegment(segment);
        }

        class UncallableDataException : Exception
        {
            public UncallableDataException()
            {
            }

            public UncallableDataException(string message)
                : base(message)
            {
            }

            public UncallableDataException(string message, Exception inner)
                : base(message, inner)
            {
            }
        }
    }
}
