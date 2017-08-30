using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.Linq;
using CanvasCommon;
using Illumina.Common.FileSystem;
using Isas.Framework.Logging;

namespace CanvasDiploidCaller
{
    class CanvasDiploidCaller
    {
        #region Members
        // Static:
        static private int MaximumCopyNumber = 10;

        // Data:
        List<CanvasSegment> Segments;
        List<SegmentPloidy> AllPloidies;
        double DiploidCoverage = 0;

        // Parameters:
        protected double MeanCoverage = 30;
        public static double CoverageWeighting = 0.6;
        private double CoverageWeightingFactor; // Computed from CoverageWeighting
        public bool IsDbsnpVcf = false;
        static protected int MinimumVariantFrequenciesForInformativeSegment = 50;
        protected int MedianHetSnpsDistance = 463; // based on NA12878 VFResults.txt.gz file
        CopyNumberOracle CNOracle = null;
        public QualityScoreParameters germlineScoreParameters;
        public int QualityFilterThreshold { get; set; } = 10;

        // File paths:
        public string TempFolder;
        CoverageModel Model;
        private readonly ILogger _logger;
        private Segments _segments;

        public CanvasDiploidCaller(ILogger logger)
        {
            _logger = logger;
        }

        #endregion

        /// <summary>
        /// Setup: Model various copy ploidies.
        /// </summary>
        public void InitializePloidies()
        {
            Console.WriteLine("{0} Initialize ploidy models...", DateTime.Now);
            this.AllPloidies = new List<SegmentPloidy>();
            CanvasCommon.Utilities.EstimateDiploidMAF(2, this.MeanCoverage);
            for (int copyNumber = 0; copyNumber <= MaximumCopyNumber; copyNumber++)
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
                        ploidy.MinorAlleleFrequency = 0.01; // should reflect sequencing error rate
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
        /// Compute the expected bin counts for each copy number, given a specified bin count for CN=2 regions
        /// </summary>
        protected static double[] GetProjectedMeanCoverage(double diploidCoverage)
        {
            double[] mu = new double[MaximumCopyNumber + 1];
            for (int count = 0; count < mu.Length; count++)
            {
                mu[count] = diploidCoverage * count / 2f;
            }
            return mu;
        }

        protected List<ModelPoint> InitializeModelPoints(CoverageModel model)
        {
            List<ModelPoint> modelPoints = new List<ModelPoint>();

            double[] mu = GetProjectedMeanCoverage(model.DiploidCoverage);
            // Refine our estimate of diploid MAF:
            //double diploidMAF = this.EstimateDiploidMAF(2, model.DiploidCoverage);

            /////////////////////////////////////////////
            // Update the parameters in each SegmentPloidy object, and construct corresponding SegmentInfo objects
            foreach (SegmentPloidy ploidy in this.AllPloidies)
            {
                ModelPoint point = new ModelPoint();
                double pureCoverage = mu[ploidy.CopyNumber];
                point.Coverage = pureCoverage;
                double pureMAF = ploidy.MinorAlleleFrequency;
                point.MAF = pureMAF;
                if (double.IsNaN(point.MAF)) point.MAF = 0;
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
        private double FitGaussians(CoverageModel model, List<SegmentInfo> segments, string debugPath = null)
        {
            List<ModelPoint> modelPoints = InitializeModelPoints(model);

            GaussianMixtureModel gmm = new GaussianMixtureModel(modelPoints, segments, this.MeanCoverage, this.CoverageWeightingFactor, 0);
            double likelihood = gmm.Fit();

            if (debugPath != null)
            {
                // write Gaussian mixture model to debugPath
                using (FileStream stream = new FileStream(debugPath, FileMode.Create, FileAccess.Write))
                using (StreamWriter writer = new StreamWriter(stream))
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

        private void AssignPloidyCallsDistance(CoverageModel model, List<SegmentInfo> segments, int medianVariantCoverage)
        {
            List<ModelPoint> modelPoints = InitializeModelPoints(model);
            foreach (CanvasSegment segment in this.Segments)
            {
                // Compute (MAF, Coverage) for this segment:
                List<double> MAF = new List<double>();
                foreach (float VF in segment.Balleles.Frequencies) MAF.Add(VF > 0.5 ? 1 - VF : VF);
                int expectedSnpDensityCutoff = (segment.End - segment.Begin) / MedianHetSnpsDistance / 2;

                double medianCoverage = CanvasCommon.Utilities.Median(segment.Counts);

                double medianMAF = -1;

                SegmentPloidy bestPloidy = null;

                if (MAF.Count >= Math.Max(10, expectedSnpDensityCutoff))
                {
                    medianMAF = Utilities.Median(MAF);
                }

                double bestDistance = double.MaxValue;
                double secondBestDistance = double.MaxValue;

                foreach (SegmentPloidy ploidy in AllPloidies)
                {
                    double diff = (ploidy.MixedCoverage - medianCoverage) * CoverageWeightingFactor;
                    double distance = diff * diff;
                    if (MAF.Count >= Math.Max(10, expectedSnpDensityCutoff))
                    {
                        diff = ploidy.MixedMinorAlleleFrequency - medianMAF;
                        distance += diff * diff;
                    }
                    if (distance < bestDistance)
                    {
                        secondBestDistance = bestDistance;
                        bestDistance = distance;
                        bestPloidy = ploidy;
                    }
                    else if (distance < secondBestDistance)
                    {
                        secondBestDistance = distance;
                    }
                }
                segment.CopyNumber = bestPloidy.CopyNumber;
                segment.ModelDistance = bestDistance;
                segment.RunnerUpModelDistance = secondBestDistance;

                segment.MajorChromosomeCount = bestPloidy.MajorChromosomeCount;
                if (MAF.Count < 10) segment.MajorChromosomeCount = null; // Don't assign MCC if we don't have variant allele frequencies
            }
        }

        /// <summary>
        /// Assign a SegmentPloidy to each CanvasSegment, based on which model matches this segment best:
        /// </summary>
        void AssignPloidyCallsGaussianMixture()
        {

            // For segments with (almost) no variants alleles at all, we'll assign them a dummy MAF, and 
            // we simply won't consider MAF when determining the closest ploidy:
            double dummyMAF = -1;

            foreach (CanvasSegment segment in this.Segments)
            {
                // Compute (MAF, Coverage) for this segment:
                List<double> MAF = new List<double>();
                foreach (float VF in segment.Balleles.Frequencies) MAF.Add(VF > 0.5 ? 1 - VF : VF);
                double medianCoverage = CanvasCommon.Utilities.Median(segment.Counts);
                double medianMAF = dummyMAF;

                SegmentPloidy bestPloidy = null;
                double bestProbability = 0;

                if (MAF.Count >= 10)
                {
                    medianMAF = Utilities.Median(MAF);
                }

                Dictionary<SegmentPloidy, double> posteriorProbabilities = GaussianMixtureModel.EMComputePosteriorProbs(AllPloidies, medianMAF, medianCoverage);
                // Find the closest ploidy. 
                foreach (SegmentPloidy ploidy in AllPloidies)
                {
                    if (bestPloidy == null || posteriorProbabilities[ploidy] > bestProbability)
                    {
                        bestProbability = posteriorProbabilities[ploidy];
                        bestPloidy = ploidy;
                    }
                }

                if (bestProbability == 0)
                {
                    // Sanity-check: If we didn't find anything with probability > 0, then fall back to the simplest possible
                    // thing: Call purely on coverage.
                    segment.CopyNumber = (int)Math.Round(2 * medianCoverage / this.DiploidCoverage);
                    segment.MajorChromosomeCount = null;
                }
                else
                {
                    segment.CopyNumber = bestPloidy.CopyNumber;
                    segment.MajorChromosomeCount = bestPloidy.MajorChromosomeCount;
                    if (MAF.Count < 10) segment.MajorChromosomeCount = null; // Don't assign MCC if we don't have variant allele frequencies
                }
            }
        }

        static public float[] AggregateCounts(ref List<CanvasSegment> segments)
        {
            List<float> diploidCounts = new List<float>();

            foreach (CanvasSegment segment in segments)
            {
                foreach (float count in segment.Counts)
                    diploidCounts.Add(count);
            }
            return diploidCounts.ToArray();
        }

        public static int AggregateVariantCoverage(ref List<CanvasSegment> segments)
        {
            var variantCoverage = segments.SelectMany(segment => segment.Balleles.TotalCoverage).ToList();
            return variantCoverage.Any() ? Utilities.Median(variantCoverage) : 0;
        }

        /// <summary>
        /// Check whether we know the CN for this segment.  Look for a known-CN interval that 
        /// covers (at least half of) this segment.  Return -1 if we don't know its CN.
        /// </summary>
        protected int GetKnownCNForSegment(CanvasSegment segment)
        {
            if (CNOracle == null) return -1;
            return CNOracle.GetKnownCNForSegment(segment);
        }

        /// <summary>
        /// Generate a table listing segments (and several features), and noting which are accurate (copy number 
        /// exactly matches truth set) or directionally accurate (copy number and truth set are both <2, both =2, or both >2)
        /// This table will become our collection of feature vectors for training q-scores!
        /// </summary>
        private void GenerateReportVersusKnownCN()
        {
            string debugPath = Path.Combine(this.TempFolder, "CallsVersusKnownCN.txt");
            using (FileStream stream = new FileStream(debugPath, FileMode.Create, FileAccess.Write))
            using (StreamWriter writer = new StreamWriter(stream))
            {
                writer.Write("#Accurate\tDirectionAccurate\t");
                writer.Write("Chr\tBegin\tEnd\tTruthSetCN\t");
                writer.Write("LogLength\tLogBinCount\tBinCount\tBinCV\tModelDistance\tRunnerUpModelDistance\t");
                writer.Write("MafCount\tMafMean\tMafCv\tLogMafCv\tCopyNumber\tMCC\t");
                writer.Write("DistanceRatio\tLogMafCount\t");
                writer.Write("ModelPurity\tModelDeviation\t");
                writer.Write("QScoreLinearFit\tQScoreGeneralizedLinearFit\tQScoreLogistic\tQScoreGermlineLogistic");
                writer.WriteLine();
                foreach (CanvasSegment segment in this.Segments)
                {
                    int CN = this.GetKnownCNForSegment(segment);
                    if (CN < 0) continue;
                    if (segment.End - segment.Begin < 5000) continue;
                    string accurateFlag = "N";
                    if (CN == segment.CopyNumber) accurateFlag = "Y";
                    string directionAccurateFlag = "N";
                    if ((CN < 2 && segment.CopyNumber < 2) ||
                        (CN == 2 && segment.CopyNumber == 2) ||
                        (CN > 2 && segment.CopyNumber > 2))
                        directionAccurateFlag = "Y";
                    writer.Write("{0}\t{1}\t", accurateFlag, directionAccurateFlag);
                    writer.Write("{0}\t{1}\t{2}\t{3}\t", segment.Chr, segment.Begin, segment.End, CN);
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
                    writer.Write("{0}\t", 100);
                    writer.Write("{0}\t", Model.Deviation);
                    double score = segment.ComputeQScore(CanvasSegment.QScoreMethod.BinCountLinearFit, germlineScoreParameters);
                    writer.Write("{0}\t", score);
                    score = segment.ComputeQScore(CanvasSegment.QScoreMethod.GeneralizedLinearFit, germlineScoreParameters);
                    writer.Write("{0}\t", score);
                    score = segment.ComputeQScore(CanvasSegment.QScoreMethod.Logistic, germlineScoreParameters);
                    writer.Write("{0}\t", score);
                    score = segment.ComputeQScore(CanvasSegment.QScoreMethod.LogisticGermline, germlineScoreParameters);
                    writer.Write("{0}\t", score);

                    writer.WriteLine();
                }
            }
            Console.WriteLine(">>> Wrote report of CNV calls versus reference calls to {0}", debugPath);
        }

        public int CallVariants(string variantFrequencyFile, string inFile, string outFile, string ploidyBedPath, string referenceFolder, string sampleName,
            string truthDataPath)
        {
            if (!string.IsNullOrEmpty(truthDataPath))
            {
                this.CNOracle = new CopyNumberOracle();
                this.CNOracle.LoadKnownCN(truthDataPath);
            }

            _segments = CanvasCommon.Segments.ReadSegments(_logger, new FileLocation(inFile));
            Segments = _segments.AllSegments.ToList();
            this.TempFolder = Path.GetDirectoryName(inFile);
            if (this.Segments.Count == 0)
            {
                Console.WriteLine("CanvasDiploidCaller: No segments loaded; no CNV calls will be made.");
                CanvasSegmentWriter.WriteSegments(outFile, this.Segments, Model?.DiploidCoverage, referenceFolder,
                    sampleName, null, null, QualityFilterThreshold, isPedigreeInfoSupplied: false);
                return 0;
            }
            PloidyInfo ploidy = null;
            if (!string.IsNullOrEmpty(ploidyBedPath)) ploidy = PloidyInfo.LoadPloidyFromBedFile(ploidyBedPath);

            // load MAF
            var allelesByChromosome = CanvasIO.ReadFrequenciesWrapper(_logger, new FileLocation(variantFrequencyFile), _segments.IntervalsByChromosome);
            _segments.AddAlleles(allelesByChromosome);
            this.MeanCoverage = allelesByChromosome.SelectMany(x => x.Value).SelectMany(y => y.TotalCoverage).Average();
            int medianVariantCoverage = AggregateVariantCoverage(ref this.Segments);

            // Create new models for different copy number states
            this.InitializePloidies();

            // Compute statistics on the copy number two regions
            float[] diploidCounts = AggregateCounts(ref this.Segments);
            DiploidCoverage = CanvasCommon.Utilities.Mean(diploidCounts);
            CoverageWeightingFactor = CoverageWeighting / DiploidCoverage;
            // new coverage model
            this.Model = new CoverageModel();
            Model.DiploidCoverage = DiploidCoverage;
            List<SegmentInfo> segments = new List<SegmentInfo>();
            foreach (CanvasSegment segment in this.Segments)
            {
                SegmentInfo info = new SegmentInfo();
                info.Segment = segment;
                List<double> MAF = new List<double>();
                foreach (float value in segment.Balleles.Frequencies) MAF.Add(value > 0.5 ? 1 - value : value);

                if (MAF.Count > 0)
                {
                    info.MAF = CanvasCommon.Utilities.Median(MAF);

                }
                else
                {
                    info.MAF = -1;
                }

                info.Coverage = CanvasCommon.Utilities.Median(segment.Counts);

                if (this.Segments.Count > 100)
                {
                    info.Weight = segment.End - segment.Begin;
                }
                else
                {
                    info.Weight = segment.BinCount;
                }
                segments.Add(info);
            }

            // Assign copy number and major chromosome count for each segment
            bool useGaussianMixtureModel = false; // For now, this is set false, since we saw weird performance on chrY (CANV-115):
            if (useGaussianMixtureModel)
            {
                // optimize model covariance
                double likelihood = FitGaussians(Model, segments);
                AssignPloidyCallsGaussianMixture();
            }
            else
            {
                AssignPloidyCallsDistance(Model, segments, medianVariantCoverage);
            }

            CanvasSegment.AssignQualityScores(this.Segments, CanvasSegment.QScoreMethod.LogisticGermline, germlineScoreParameters);

            // Merge neighboring segments that got the same copy number call.
            // merging segments requires quality scores so we do it after quality scores have been assigned
            var mergedSegments = CanvasSegment.MergeSegments(Segments);
            // recalculating qscores after merging segments improves performance!
            CanvasSegment.AssignQualityScores(mergedSegments, CanvasSegment.QScoreMethod.LogisticGermline, germlineScoreParameters);
            CanvasSegment.FilterSegments(QualityFilterThreshold, mergedSegments); 

            List<string> extraHeaders = new List<string>();
            string coverageOutputPath = SingleSampleCallset.GetCoverageAndVariantFrequencyOutputPath(outFile);
            CanvasSegment.WriteCoveragePlotData(mergedSegments, Model.DiploidCoverage, ploidy, coverageOutputPath, referenceFolder);

            if (this.CNOracle != null)
            {
                this.GenerateReportVersusKnownCN();
            }

            if (ploidy != null && !string.IsNullOrEmpty(ploidy.HeaderLine)) extraHeaders.Add(ploidy.HeaderLine);

            CanvasSegmentWriter.WriteSegments(outFile, mergedSegments, Model.DiploidCoverage, referenceFolder, sampleName,
                extraHeaders, ploidy, QualityFilterThreshold, isPedigreeInfoSupplied: false);
            return 0;
        }
    }
}
