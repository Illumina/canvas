using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using CanvasCommon;
using Illumina.Common.FileSystem;
using Isas.Framework.Logging;

namespace CanvasDiploidCaller
{
    internal class CanvasDiploidCaller
    {
        #region Members
        // Static:
        private const int MaximumCopyNumber = 10;

        // Data:
        private List<CanvasSegment> _allSegments;

        private List<SegmentPloidy> _allPloidies;
        private double _diploidCoverage;

        // Parameters:
        protected double MeanCoverage = 30;
        public static double CoverageWeighting = 0.6;
        private double _coverageWeightingFactor; // Computed from CoverageWeighting
        protected static int MinimumVariantFrequenciesForInformativeSegment = 50;
        protected int MedianHetSnpsDistance = 463; // based on NA12878 VFResults.txt.gz file
        CopyNumberOracle _cnOracle;
        private readonly QualityScoreParameters _germlineScoreParameters;
        public int QualityFilterThreshold { get; set; } = 10;

        // File paths:
        public string TempFolder;

        private CoverageModel _model;
        private readonly ILogger _logger;
        private Segments _segments;

        public CanvasDiploidCaller(ILogger logger, QualityScoreParameters qscoreParametersJson)
        {
            _logger = logger;
            _germlineScoreParameters = qscoreParametersJson;
        }

        #endregion

        /// <summary>
        /// Setup: Model various copy ploidies.
        /// </summary>
        public void InitializePloidies()
        {
            Console.WriteLine("{0} Initialize ploidy models...", DateTime.Now);
            _allPloidies = new List<SegmentPloidy>();
            Utilities.EstimateDiploidMAF(2, MeanCoverage);
            for (int copyNumber = 0; copyNumber <= MaximumCopyNumber; copyNumber++)
            {
                for (int majorCount = copyNumber; majorCount * 2 >= copyNumber; majorCount--)
                {
                    SegmentPloidy ploidy = new SegmentPloidy
                    {
                        CopyNumber = copyNumber,
                        MajorChromosomeCount = majorCount,
                        Index = _allPloidies.Count
                    };
                    _allPloidies.Add(ploidy);
                    if (copyNumber == 0)
                    {
                        ploidy.MinorAlleleFrequency = 0.01; // should reflect sequencing error rate
                        continue;
                    }
                    float variantFrequency = majorCount / (float)copyNumber;
                    ploidy.MinorAlleleFrequency = variantFrequency < 0.5 ? variantFrequency : 1 - variantFrequency;
                    if (majorCount * 2 == copyNumber)
                    {
                        ploidy.MinorAlleleFrequency = Utilities.EstimateDiploidMAF(copyNumber, MeanCoverage);
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
            foreach (SegmentPloidy ploidy in _allPloidies)
            {
                ModelPoint point = new ModelPoint();
                double pureCoverage = mu[ploidy.CopyNumber];
                point.Coverage = pureCoverage;
                double pureMaf = ploidy.MinorAlleleFrequency;
                point.Maf = pureMaf;
                if (double.IsNaN(point.Maf)) point.Maf = 0;
                point.Ploidy = ploidy;
                modelPoints.Add(point);
                point.CopyNumber = ploidy.CopyNumber;
                ploidy.MixedMinorAlleleFrequency = point.Maf;
                ploidy.MixedCoverage = point.Coverage;
            }

            return modelPoints;
        }

        private void AssignPloidyCallsDistance(CoverageModel model)
        {
            InitializeModelPoints(model);
            foreach (CanvasSegment segment in _allSegments)
            {
                // Compute (MAF, Coverage) for this segment:
                List<double> mafs = new List<double>();
                foreach (float variantFrequency in segment.Balleles.Frequencies) mafs.Add(variantFrequency > 0.5 ? 1 - variantFrequency : variantFrequency);
                int expectedSnpDensityCutoff = (segment.Length) / MedianHetSnpsDistance / 2;

                double medianCoverage = Utilities.Median(segment.Counts);

                double medianMaf = -1;

                SegmentPloidy bestPloidy = null;

                if (mafs.Count >= Math.Max(10, expectedSnpDensityCutoff))
                {
                    medianMaf = Utilities.Median(mafs);
                }

                double bestDistance = double.MaxValue;
                double secondBestDistance = double.MaxValue;

                foreach (SegmentPloidy ploidy in _allPloidies)
                {
                    double diff = (ploidy.MixedCoverage - medianCoverage) * _coverageWeightingFactor;
                    double distance = diff * diff;
                    if (mafs.Count >= Math.Max(10, expectedSnpDensityCutoff))
                    {
                        diff = ploidy.MixedMinorAlleleFrequency - medianMaf;
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
                if (bestPloidy != null)
                {
                    segment.CopyNumber = bestPloidy.CopyNumber;
                    segment.MajorChromosomeCount = bestPloidy.MajorChromosomeCount;
                }
                segment.ModelDistance = bestDistance;
                segment.RunnerUpModelDistance = secondBestDistance;

                if (mafs.Count < 10) segment.MajorChromosomeCount = null; // Don't assign MCC if we don't have variant allele frequencies
            }
        }

        public static float[] AggregateCounts(ref List<CanvasSegment> segments)
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
        protected int GetKnownCopyNumberForSegment(CanvasSegment segment)
        {
            if (_cnOracle == null) return -1;
            return _cnOracle.GetKnownCNForSegment(segment);
        }

        /// <summary>
        /// Generate a table listing segments (and several features), and noting which are accurate (copy number 
        /// exactly matches truth set) or directionally accurate (copy number and truth set are both &lt;2, both =2, or both &gt;2)
        /// This table will become our collection of feature vectors for training q-scores!
        /// </summary>
        private void GenerateReportVersusKnownCopyNumber()
        {
            string debugPath = Path.Combine(TempFolder, "CallsVersusKnownCN.txt");
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
                foreach (CanvasSegment segment in _allSegments)
                {
                    int copyNumber = GetKnownCopyNumberForSegment(segment);
                    if (copyNumber < 0) continue;
                    if (segment.Length < 5000) continue;
                    string accurateFlag = "N";
                    if (copyNumber == segment.CopyNumber) accurateFlag = "Y";
                    string directionAccurateFlag = "N";
                    if ((copyNumber < 2 && segment.CopyNumber < 2) ||
                        (copyNumber == 2 && segment.CopyNumber == 2) ||
                        (copyNumber > 2 && segment.CopyNumber > 2))
                        directionAccurateFlag = "Y";
                    writer.Write("{0}\t{1}\t", accurateFlag, directionAccurateFlag);
                    writer.Write("{0}\t{1}\t{2}\t{3}\t", segment.Chr, segment.Begin, segment.End, copyNumber);
                    writer.Write("{0}\t", Math.Log(segment.Length));
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
                    writer.Write("{0}\t", _model.Deviation);
                    double score = segment.ComputeQScore(CanvasSegment.QScoreMethod.BinCountLinearFit, _germlineScoreParameters);
                    writer.Write("{0}\t", score);
                    score = segment.ComputeQScore(CanvasSegment.QScoreMethod.GeneralizedLinearFit, _germlineScoreParameters);
                    writer.Write("{0}\t", score);
                    score = segment.ComputeQScore(CanvasSegment.QScoreMethod.Logistic, _germlineScoreParameters);
                    writer.Write("{0}\t", score);
                    score = segment.ComputeQScore(CanvasSegment.QScoreMethod.LogisticGermline, _germlineScoreParameters);
                    writer.Write("{0}\t", score);

                    writer.WriteLine();
                }
            }
            Console.WriteLine(">>> Wrote report of CNV calls versus reference calls to {0}", debugPath);
        }

        public int CallVariants(string variantFrequencyFile, string inFile, string outFile, string ploidyVcfPath, string referenceFolder, string sampleName,
            string truthDataPath)
        {
            if (!string.IsNullOrEmpty(truthDataPath))
            {
                _cnOracle = new CopyNumberOracle();
                _cnOracle.LoadKnownCN(truthDataPath);
            }

            _segments = Segments.ReadSegments(_logger, new FileLocation(inFile));
            _allSegments = _segments.AllSegments.ToList();
            TempFolder = Path.GetDirectoryName(inFile);
            if (_allSegments.Count == 0)
            {
                Console.WriteLine("CanvasDiploidCaller: No segments loaded; no CNV calls will be made.");
                CanvasSegmentWriter.WriteSegments(outFile, _allSegments, _model?.DiploidCoverage, referenceFolder,
                    sampleName, null, null, QualityFilterThreshold, isPedigreeInfoSupplied: false);
                return 0;
            }
            PloidyInfo ploidy = null;
            if (!string.IsNullOrEmpty(ploidyVcfPath)) ploidy = PloidyInfo.LoadPloidyFromVcfFileNoSampleId(ploidyVcfPath);

            // load MAF
            var allelesByChromosome = CanvasIO.ReadFrequenciesWrapper(_logger, new FileLocation(variantFrequencyFile), _segments.IntervalsByChromosome);
            _segments.AddAlleles(allelesByChromosome);
            MeanCoverage = allelesByChromosome.SelectMany(x => x.Value).SelectMany(y => y.TotalCoverage).Average();
            AggregateVariantCoverage(ref _allSegments);

            // Create new models for different copy number states
            InitializePloidies();

            // Compute statistics on the copy number two regions
            float[] diploidCounts = AggregateCounts(ref _allSegments);
            _diploidCoverage = Utilities.Mean(diploidCounts);
            _coverageWeightingFactor = CoverageWeighting / _diploidCoverage;
            // new coverage model
            _model = new CoverageModel { DiploidCoverage = _diploidCoverage };
            List<SegmentInfo> segments = new List<SegmentInfo>();
            foreach (CanvasSegment segment in _allSegments)
            {
                SegmentInfo info = new SegmentInfo { Segment = segment };
                List<double> mafs = new List<double>();
                foreach (float value in segment.Balleles.Frequencies) mafs.Add(value > 0.5 ? 1 - value : value);

                if (mafs.Count > 0)
                {
                    info.Maf = Utilities.Median(mafs);

                }
                else
                {
                    info.Maf = -1;
                }

                info.Coverage = Utilities.Median(segment.Counts);

                info.Weight = _allSegments.Count > 100 ? segment.Length : segment.BinCount;
                segments.Add(info);
            }

            AssignPloidyCallsDistance(_model);

            CanvasSegment.AssignQualityScores(_allSegments, CanvasSegment.QScoreMethod.LogisticGermline, _germlineScoreParameters);

            // Merge neighboring segments that got the same copy number call.
            // merging segments requires quality scores so we do it after quality scores have been assigned
            var mergedSegments = CanvasSegment.MergeSegments(_allSegments);
            // recalculating qscores after merging segments improves performance!

            CanvasSegment.AssignQualityScores(mergedSegments, CanvasSegment.QScoreMethod.LogisticGermline, _germlineScoreParameters);
            CanvasSegment.SetFilterForSegments(QualityFilterThreshold, mergedSegments, CanvasFilter.SegmentSizeCutoff); 

            List<string> extraHeaders = new List<string>();
            var coverageOutputPath = SingleSampleCallset.GetCoverageAndVariantFrequencyOutputPath(outFile);
            CanvasSegment.WriteCoveragePlotData(mergedSegments, _model.DiploidCoverage, ploidy, coverageOutputPath, referenceFolder);

            if (_cnOracle != null)
            {
                GenerateReportVersusKnownCopyNumber();
            }

            if (!string.IsNullOrEmpty(ploidy?.HeaderLine)) extraHeaders.Add(ploidy.HeaderLine);

            CanvasSegmentWriter.WriteSegments(outFile, mergedSegments, _model.DiploidCoverage, referenceFolder, sampleName,
                extraHeaders, ploidy, QualityFilterThreshold, isPedigreeInfoSupplied: false);
            return 0;
        }
    }
}
