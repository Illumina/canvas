using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NDesk.Options;
using System.IO;
using SequencingFiles;
using CanvasCommon;

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
        protected float MeanCoverage = 30;
        public static double CoverageWeighting = 0.4;
        private double CoverageWeightingFactor; // Computed from CoverageWeighting
        public bool IsDbsnpVcf = false;
        static protected int MinimumVariantFrequenciesForInformativeSegment = 50;

        // File paths:
        public string TempFolder;

        #endregion

        /// <summary>
        /// Setup: Model various copy ploidies.
        /// </summary>
        public void InitializePloidies()
        {
            Console.WriteLine("{0} Initialize ploidy models...", DateTime.Now);
            this.AllPloidies = new List<SegmentPloidy>();
            double diploidPredictedMAF = CanvasCommon.Utilities.EstimateDiploidMAF(2, this.MeanCoverage);
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
        /// Compute the expected bin counts for each copy number, given a specified bin count for CN=2 regions, and assuming
        /// pure tumor data.  
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
            double diploidMAF = this.AllPloidies[3].MinorAlleleFrequency; /// %%% Magic number!
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
                point.MAF =  pureMAF;
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

            GaussianMixtureModel gmm = new GaussianMixtureModel(modelPoints, segments, this.MeanCoverage, this.CoverageWeightingFactor);
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

                // Sanity-check: If we didn't find anything with probability > 0, then fall back to the simplest possible
                // thing: Call purely on coverage.
                if (bestProbability == 0)
                {
                    segment.CopyNumber = (int)Math.Round(2 * medianCoverage / this.DiploidCoverage);
                    segment.MajorChromosomeCount = (int)Math.Ceiling(segment.CopyNumber / 2f);
                }
                else
                {
                    segment.CopyNumber = bestPloidy.CopyNumber;
                    segment.PloidyIndex = bestPloidy.ID;
                    segment.MajorChromosomeCount = bestPloidy.MajorChromosomeCount;
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


        public int CallVariants(string variantFrequencyFile, string inFile, string outFile, string ploidyBedPath, string referenceFolder, string sampleName)
        {
            this.Segments = CanvasSegment.ReadSegments(inFile);
            this.TempFolder = Path.GetDirectoryName(inFile);
            if (this.Segments.Count == 0)
            {
                Console.WriteLine("CanvasDiploidCaller: No segments loaded; no CNV calls will be made.");
                CanvasSegment.WriteSegments(outFile, this.Segments, referenceFolder, sampleName, null, false, null, false);
                return 0;
            }
            PloidyInfo ploidy = null;
            if (!string.IsNullOrEmpty(ploidyBedPath)) ploidy = PloidyInfo.LoadPloidyFromBedFile(ploidyBedPath);
            double diploidCount = CanvasSegment.ExpectedCount(this.Segments);

            // load MAF
            this.MeanCoverage = CanvasIO.LoadVariantFrequencies(variantFrequencyFile, this.Segments);

            // Create new models for different copy number states
            this.InitializePloidies();

            // Compute statistics on the copy number two regions
            float[] diploidCounts = AggregateCounts(ref this.Segments);
            DiploidCoverage = CanvasCommon.Utilities.Mean(diploidCounts);
            CoverageWeightingFactor = CoverageWeighting / DiploidCoverage;

            // new coverage model
            CoverageModel model = new CoverageModel();
            model.DiploidCoverage = DiploidCoverage;
            List<SegmentInfo> segments = new List<SegmentInfo>();
            foreach (CanvasSegment segment in this.Segments)
            {
                SegmentInfo info = new SegmentInfo();
                info.Segment = segment;
                List<double> MAF = new List<double>();
                foreach (float value in segment.VariantFrequencies) MAF.Add(value > 0.5 ? 1 - value : value);
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
            // optimize model covariance
            double likelihood = FitGaussians(model, segments);
            // Assign copy number and major chromosome count for each segment
            AssignPloidyCallsGaussianMixture();

            // Merge neighboring segments that got the same copy number call.
            CanvasSegment.MergeSegments(ref this.Segments);
            CanvasSegment.AssignQualityScores(this.Segments, CanvasSegment.QScoreMethod.BinCountLinearFit);
            List<string> extraHeaders = new List<string>();

            string coverageOutputPath = CanvasCommon.Utilities.GetCoverageAndVariantFrequencyOutputPath(outFile);
            CanvasSegment.WriteCoveragePlotData(this.Segments, model.DiploidCoverage, ploidy, coverageOutputPath, referenceFolder);

            if (ploidy != null && !string.IsNullOrEmpty(ploidy.HeaderLine)) extraHeaders.Add(ploidy.HeaderLine);
            CanvasSegment.WriteSegments(outFile, this.Segments, referenceFolder, sampleName, extraHeaders, true, ploidy, true, false);


            return 0;
        }
    }
}
