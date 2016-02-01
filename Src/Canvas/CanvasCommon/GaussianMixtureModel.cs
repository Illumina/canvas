using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CanvasCommon
{
    
    public class GaussianMixtureModel
    {
        #region Members
        List<ModelPoint> ModelPoints;
        List<SegmentInfo> Segments;
        private double MeanCoverage;
        private double CoverageWeightingFactor;
        private double KnearestNeighbourCutoff;

        // parameters
        public bool uncorrelated = true; // Assume MAF and Coverage are independent/uncorrelated in FFPEMode (always false for now)
        private const double EMPosteriorProbThres = 0.01; // Controls whether a segment contributes to Mu and Sigma estimates
        private const double EMOmegaThres = 0.01; // Controls when to update means
        private const double EMLikelihoodThres = 1; // Controls when to update means
        #endregion

        public GaussianMixtureModel(List<ModelPoint> modelPoints, List<SegmentInfo> segments,
            double meanCoverage, double coverageWeightingFactor, double knearestNeighbourCutoff) 
        {
            ModelPoints = modelPoints;
            Segments = segments;
            MeanCoverage = meanCoverage;
            CoverageWeightingFactor = coverageWeightingFactor;
            KnearestNeighbourCutoff = knearestNeighbourCutoff;
        }

        public double Fit() 
        {
            double likelihood = FitGaussians(this.ModelPoints, this.Segments);
            EMComputeGaussianMeans(this.ModelPoints, this.Segments);
            foreach (var modelPoint in this.ModelPoints) { modelPoint.Ploidy.Sigma = null; }
            likelihood = FitGaussians(this.ModelPoints, this.Segments);
            return likelihood;
        }

        private void EMInitializeParameters(List<ModelPoint> modelPoints)
        {
            const int Dim = 2;

            foreach (ModelPoint modelPoint in modelPoints)
            {
                modelPoint.Ploidy.Omega = 1.0 / modelPoints.Count;
                modelPoint.Ploidy.Mu = new double[Dim];
                modelPoint.Ploidy.Mu[0] = modelPoint.MAF;
                modelPoint.Ploidy.Mu[1] = modelPoint.Coverage;

                // initialize covariance matrix
                modelPoint.Ploidy.Sigma = new double[Dim][];
                for (int j = 0; j < Dim; j++) { modelPoint.Ploidy.Sigma[j] = new double[Dim]; }
                modelPoint.Ploidy.Sigma[0][0] = 0.01;
                modelPoint.Ploidy.Sigma[0][1] = 0.0;
                modelPoint.Ploidy.Sigma[1][0] = 0.0;
                modelPoint.Ploidy.Sigma[1][1] = 0.01 / (this.CoverageWeightingFactor * this.CoverageWeightingFactor);
            }
        }

        private static double Sigma(double intensity, double mu, double sigma2)
        {
            double likelihood = 0;

            double diff = intensity - mu;
            double tempval = -0.5 * (diff * diff / sigma2);

            if (!Double.IsNaN(tempval)) //check for nans
            {
                likelihood = 1.0 / (Math.Sqrt(2.0 * Math.PI * sigma2)) * Math.Exp(tempval);
            }
            if (Double.IsNaN(likelihood)) { likelihood = 0; }

            return likelihood;
        }

        /// <summary>
        /// Calculate the Sigma function for use in the EM calculation. Assumes 2-D data.
        /// </summary>
        /// <param name="intensityX"></param>
        /// <param name="intensityY"></param>
        /// <param name="mu">means</param>
        /// <param name="sigma">covariance matrix</param>
        /// <returns></returns>
        private static double Sigma(double intensityX, double intensityY, double[] mu, double[][] sigma)
        {
            if (intensityX == -1 ) // dummy MAF
            {
                Sigma(intensityY, mu[1], sigma[1][1]);
            }

            double likelihood = 0;
            //calculate the determinant of Sigma
            double det = sigma[0][0] * sigma[1][1] - sigma[0][1] * sigma[1][0];

            //now we need mu*invSigma*mu'
            double diffX = intensityX - mu[0];
            double diffY = intensityY - mu[1];
            double diffProduct = diffX * diffY;
            double tempval = -0.5 * (sigma[1][1] / det * diffX * diffX -
                sigma[0][1] / det * diffProduct -
                sigma[1][0] / det * diffProduct +
                sigma[0][0] / det * diffY * diffY);

            if (!Double.IsNaN(tempval)) //check for nans
            {
                likelihood = 1.0 / (2.0 * Math.PI * Math.Sqrt(det)) * Math.Exp(tempval);
            }
            if (Double.IsNaN(likelihood)) { likelihood = 0; }

            return likelihood;
        }

        private void EMComputePosteriorProbs(List<ModelPoint> modelPoints, SegmentInfo segment)
        {
            double tempsum1 = 0;
            Dictionary<ModelPoint, double> temp = new Dictionary<ModelPoint, double>();
            foreach (var modelPoint in modelPoints)
            {
                temp[modelPoint] = modelPoint.Ploidy.Omega * Sigma(segment.MAF, segment.Coverage, modelPoint.Ploidy.Mu, modelPoint.Ploidy.Sigma);
                tempsum1 += temp[modelPoint];
            }

            if (segment.PosteriorProbs == null) { segment.PosteriorProbs = new Dictionary<ModelPoint, double>(); }
            int bestCluster = 0;
            double bestProb = 0;
            foreach (var modelPoint in modelPoints)
            {
                segment.PosteriorProbs[modelPoint] = temp[modelPoint] / tempsum1;
                if (segment.PosteriorProbs[modelPoint] > bestProb) { 
                    bestCluster = modelPoint.Cluster.Value;
                    bestProb = temp[modelPoint] / tempsum1;
                }
                if (Double.IsNaN(segment.PosteriorProbs[modelPoint])) { segment.PosteriorProbs[modelPoint] = 0; }
            }
            segment.Cluster = bestCluster;
        }

        /// <summary>
        /// For each segment, compute the posterior probability of the segment belonging to each model point
        /// </summary>
        /// <param name="modelPoints"></param>
        /// <param name="segments"></param>
        private void EMComputePosteriorProbs(List<ModelPoint> modelPoints, List<SegmentInfo> segments)
        {
            foreach (var segment in segments)
            {
                if (segment.KnearestNeighbour > this.KnearestNeighbourCutoff) segment.Cluster = -1;
                else EMComputePosteriorProbs(modelPoints, segment);
            }
        }

        public static Dictionary<CanvasCommon.SegmentPloidy, double> EMComputePosteriorProbs(
            IEnumerable<CanvasCommon.SegmentPloidy> ploidies, double maf, double coverage)
        {
            double tempsum1 = 0;
            Dictionary<CanvasCommon.SegmentPloidy, double> temp = new Dictionary<CanvasCommon.SegmentPloidy, double>();
            foreach (var ploidy in ploidies)
            {
                temp[ploidy] = ploidy.Omega * Sigma(maf, coverage, ploidy.Mu, ploidy.Sigma);
                tempsum1 += temp[ploidy];
            }

            Dictionary<CanvasCommon.SegmentPloidy, double> PosteriorProbs = new Dictionary<CanvasCommon.SegmentPloidy, double>();
            foreach (var ploidy in ploidies)
            {
                PosteriorProbs[ploidy] = temp[ploidy] / tempsum1;

                if (Double.IsNaN(PosteriorProbs[ploidy])) { PosteriorProbs[ploidy] = 0; }
            }

            return PosteriorProbs;
        }

        private void EMComputeOmegas(List<ModelPoint> modelPoints, List<SegmentInfo> segments)
        {
            foreach (var modelPoint in modelPoints)
            {
                double sumWeights = 0;
                modelPoint.Ploidy.Omega = 0;
                foreach (var segment in segments)
                {
                    if (segment.Cluster != -1)
                    // weight segments by segment.Weight
                    {
                        sumWeights += segment.Weight;
                        modelPoint.Ploidy.Omega += segment.PosteriorProbs[modelPoint] * segment.Weight;
                    }
                }
                modelPoint.Ploidy.Omega /= sumWeights;
            }
        }

        /// <summary>
        /// Compute eigen values of a 2-D symmetric matrix
        /// </summary>
        /// <param name="symMat">2 x 2 symmetric matrix</param>
        /// <returns>eigen values</returns>
        private double[] EigenValues2D(double[][] symMat)
        {
            double[] eigenValues = new double[2];
            double b = -(symMat[0][0] + symMat[1][1]);
            double c = symMat[0][0] * symMat[1][1] - symMat[0][1] * symMat[0][1];
            double d = b * b - 4 * c;
            if (d >= 0)
            {
                d = Math.Sqrt(d);
            }
            else // not possible theoretically
            {
                d = 0;
            }

            eigenValues[0] = (-b + d) / 2;
            eigenValues[1] = (-b - d) / 2;

            return eigenValues;
        }

        private bool IsPositiveSemiDefinite(double[][] symMat)
        {
            double[] eigenValues = EigenValues2D(symMat);

            for (int i = 0; i < 2; i++)
            {
                if (eigenValues[i] < 0) { return false; }
            }

            return true;
        }

        private void Scale2DMatrix(double[][] mat, double scale)
        {
            for (int i = 0; i < mat.Length; i++)
            {
                for (int j = 0; j < mat[i].Length; j++)
                {
                    mat[i][j] *= scale;
                }
            }
        }

        /// <summary>
        /// make sure a component doesn't "invade" the other components
        /// </summary>
        /// <param name="modelPoints"></param>
        private void EMScaleCovariancesPairwise(List<ModelPoint> modelPoints)
        {
            foreach (var m1 in modelPoints)
            {
                ModelPoint maxM = null;
                double maxProb = 0;
                foreach (var m2 in modelPoints)
                {
                    if (m2 == m1) { continue; }
                    // pretend m1 is a segment
                    double prob = m2.Ploidy.Omega * Sigma(m1.MAF, m1.Coverage, m2.Ploidy.Mu, m2.Ploidy.Sigma);
                    if (prob > maxProb)
                    {
                        maxProb = prob;
                        maxM = m2;
                    }
                }
                if (maxProb > 0)
                {
                    double[][] s1 = m1.Ploidy.Sigma;
                    double det1 = s1[0][0] * s1[1][1] - s1[0][1] * s1[1][0];
                    double[][] s2 = maxM.Ploidy.Sigma;
                    double det2 = s2[0][0] * s2[1][1] - s2[0][1] * s2[1][0];
                    if (det1 <= 1E-7 || det2 <= 1E-7) { continue; }

                    double ratio = det1 > det2 ? det1 / det2 : det2 / det1;
                    if (ratio < 4) { continue; }
                    if (det1 > det2)
                    {
                        Scale2DMatrix(s1, 0.8);
                        Scale2DMatrix(s2, 1.1);
                    }
                    else
                    {
                        Scale2DMatrix(s2, 0.8);
                        Scale2DMatrix(s1, 1.1);
                    }
                }
            }
        }

        private void EMComputeGaussianCovariances(List<ModelPoint> modelPoints, List<SegmentInfo> segments,
            double postProbThres = EMPosteriorProbThres)
        {
            int Dim = modelPoints.First().Ploidy.Mu.Length;

            double[] tempvec = new double[Dim];
            double[][] sumMat = new double[Dim][];
            for (int ix = 0; ix < Dim; ix++)
            {
                sumMat[ix] = new double[Dim];
            }

            foreach (var modelPoint in modelPoints)
            {
                for (int ix = 0; ix < Dim; ix++)
                {
                    for (int iy = 0; iy < Dim; iy++)
                    {
                        sumMat[ix][iy] = 0;
                    }
                }

                double sumWeights = 0;
                foreach (var segment in segments)
                {
                    // weight segments by segment.Weight
                    // each segment represents segment.Weight points
                    if (segment.Cluster == -1) { continue; }
                    if (segment.PosteriorProbs[modelPoint] < postProbThres) { continue; }
                    double weight = segment.PosteriorProbs[modelPoint] * segment.Weight;
                    sumWeights += weight;

                    tempvec[0] = segment.MAF - modelPoint.Ploidy.Mu[0];
                    tempvec[1] = segment.Coverage - modelPoint.Ploidy.Mu[1];

                    sumMat[0][0] += weight * tempvec[0] * tempvec[0];
                    sumMat[0][1] += weight * tempvec[0] * tempvec[1];
                    sumMat[1][0] += weight * tempvec[1] * tempvec[0];
                    sumMat[1][1] += weight * tempvec[1] * tempvec[1];
                }
                if (sumWeights > 0)
                {
                    modelPoint.Ploidy.Sigma[0][0] = sumMat[0][0] / sumWeights;
                    if (uncorrelated) // assume uncorrelated
                    {
                        modelPoint.Ploidy.Sigma[0][1] = 0;
                        modelPoint.Ploidy.Sigma[1][0] = 0;
                    }
                    else
                    {
                        modelPoint.Ploidy.Sigma[0][1] = sumMat[0][1] / sumWeights;
                        modelPoint.Ploidy.Sigma[1][0] = sumMat[1][0] / sumWeights;
                    }
                    modelPoint.Ploidy.Sigma[1][1] = sumMat[1][1] / sumWeights;
                }

                if (modelPoint.Ploidy.Sigma[0][0] < 1E-7) { modelPoint.Ploidy.Sigma[0][0] = 1E-7; }
                if (modelPoint.Ploidy.Sigma[1][1] < 1E-7) { modelPoint.Ploidy.Sigma[1][1] = 1E-7; }


                if (!IsPositiveSemiDefinite(modelPoint.Ploidy.Sigma)) // make sure the covariance matrix is positive-semidefinite
                {
                    modelPoint.Ploidy.Sigma[0][0] = 0.01;
                    modelPoint.Ploidy.Sigma[0][1] = 0.0;
                    modelPoint.Ploidy.Sigma[1][0] = 0.0;
                    modelPoint.Ploidy.Sigma[1][1] = 0.01 / (this.CoverageWeightingFactor * this.CoverageWeightingFactor);
                }
            }

            // make sure a component doesn't "invade" the other components
            EMScaleCovariancesPairwise(modelPoints);
        }

        private void EMComputeGaussianMeans(List<ModelPoint> modelPoints, List<SegmentInfo> segments,
            double postProbThres = EMPosteriorProbThres, double omegaThres = EMOmegaThres,
            double likelihoodThres = EMLikelihoodThres)
        {
            foreach (var modelPoint in modelPoints)
            {
               if (modelPoint.Ploidy.Omega < omegaThres) { continue; }

                double[] tempsum = new double[modelPoint.Ploidy.Mu.Length];
                for (int i = 0; i < tempsum.Length; i++) { tempsum[i] = 0; }

                double sumWeights = 0;
                foreach (var segment in segments)
                {
                    // weight segments by segment.Weight
                    // each segment represents segment.Weight points
                    if (segment.Cluster == -1) { continue; }
                    if (segment.PosteriorProbs[modelPoint] < postProbThres) { continue; }
                    double weight = segment.PosteriorProbs[modelPoint] * segment.Weight;
                    sumWeights += weight;
                    tempsum[0] += weight * segment.MAF;
                    tempsum[1] += weight * segment.Coverage;
                }
                double mu0 = tempsum[0] / sumWeights;
                double mu1 = tempsum[1] / sumWeights;

                modelPoint.Ploidy.Mu[0] = mu0;
                modelPoint.Ploidy.Mu[1] = mu1;

            }
        }

        private double EMComputeLikelihood(List<ModelPoint> modelPoints, List<SegmentInfo> segments)
        {
            double likelihood = 0;

            foreach (var segment in segments)
            {
                if (segment.Cluster == -1) { continue; }
                double temp = 0;
                foreach (var modelPoint in modelPoints)
                {
                    // do not consider segments with no MAFs during likelyhood calculation
                    if (segment.MAF == -1)
                        temp += modelPoint.Ploidy.Omega;
                    else
                        temp += modelPoint.Ploidy.Omega * Sigma(segment.MAF, segment.Coverage, modelPoint.Ploidy.Mu, modelPoint.Ploidy.Sigma);
                }
                likelihood += Math.Log(temp) * segment.Weight; // each segment represents segment.Weight points
            }

            return likelihood / segments.Select(s => s.Weight).Sum();
        }

        /// <summary>
        /// http://en.wikipedia.org/wiki/Bayesian_information_criterion
        /// </summary>
        /// <param name="modelPoints"></param>
        /// <param name="segments"></param>
        /// <param name="omegaThres"></param>
        /// <returns></returns>
        private double EMComputeBIC(List<ModelPoint> modelPoints, List<SegmentInfo> segments,
            double omegaThres = 1E-10, double sigmaThres = 1E-7)
        {
            double likelihood = 0; // log likelihood

            foreach (var segment in segments)
            {
                if (segment.Cluster == -1) { continue; }
                double temp = 0;
                foreach (var modelPoint in modelPoints)
                {
                    if (modelPoint.Ploidy.Omega <= omegaThres) { continue; }
                    temp += modelPoint.Ploidy.Omega * Sigma(segment.MAF, segment.Coverage, modelPoint.Ploidy.Mu, modelPoint.Ploidy.Sigma);
                }
                likelihood += Math.Log(temp) * segment.Weight; // each segment represents segment.Weight points
            }

            int k = 0;
            foreach (var modelPoint in modelPoints) 
            {
                if (modelPoint.Ploidy.Omega <= omegaThres) { continue; }
                k += 6; // number of parameters: omega 1, mu 2, sigma 3
                if (modelPoint.Ploidy.Sigma[0][0] <= sigmaThres || modelPoint.Ploidy.Sigma[1][1] <= sigmaThres) // omega 1, mu 1, sigma 1
                {   // modelPoint.Ploidy.Sigma[0][0] <= sigmaThres && modelPoint.Ploidy.Sigma[1][1] <= sigmaThres: omega 1, mu 2
                    k -= 3;
                }
                else if (Math.Abs(modelPoint.Ploidy.Sigma[0][1]) <= sigmaThres) 
                { // omega 1, mu 2, sigma 2
                    k -= 1;
                }
            }

            double n = segments.Select(s => s.Weight).Sum(); // number of data points
            double bic = -2 * likelihood + k * (Math.Log(n) - Math.Log(2 * Math.PI));

            return bic / n;
        }


        private double FitGaussians(List<ModelPoint> modelPoints, List<SegmentInfo> segments)
        {
            const int numIterations = 20; //don't go past this many optimizations
            const double likelihoodCutoff = 0.000025; //if likelihood doesn't change by this much, terminate optimization

            EMInitializeParameters(modelPoints);

            double likelihood = 0;
            double prevLikelihood = -1;

            try
            {
                // enter optimization loop
                for (int iterationCount = 0; iterationCount < numIterations; iterationCount++)
                {
                    EMComputePosteriorProbs(modelPoints, segments);

                    // update model parameters //
                    EMComputeOmegas(modelPoints, segments);

                    EMComputeGaussianCovariances(modelPoints, segments); 
                    /////////////////////////////

                    // calculate overall likelihood
                    likelihood = EMComputeLikelihood(modelPoints, segments);
                    if (Math.Abs(likelihood - prevLikelihood) < likelihoodCutoff && iterationCount > 1) { break; }
                    prevLikelihood = likelihood;
                }
            }
            catch (Exception ex)
            {
                Console.Error.WriteLine("Error fitting Gaussian mixture model: {0}", ex.Message);
                foreach (var modelPoint in modelPoints)
                {
                    modelPoint.Ploidy.Mu = null;
                    modelPoint.Ploidy.Sigma = null;
                }
            }

            return likelihood;
        }

        /// <summary>
        /// Run Gaussian expectation–maximization algorithm for segment clustering
        /// </summary>
        public double runExpectationMaximization()
        {
            const int numIterations = 30; //don't go past this many optimizations
            const double likelihoodCutoff = 0.000025; //if likelihood doesn't change by this much, terminate optimization

            EMInitializeParameters(this.ModelPoints);

            double likelihood = 0;
            double prevLikelihood = -1;

            try
            {
                // enter optimization loop
                for (int iterationCount = 0; iterationCount < numIterations; iterationCount++)
                {
                    EMComputePosteriorProbs(this.ModelPoints, this.Segments);

                    // update model parameters 
                    EMComputeOmegas(this.ModelPoints, this.Segments);
                    EMComputeGaussianMeans(this.ModelPoints, this.Segments);
                    EMComputeGaussianCovariances(this.ModelPoints, this.Segments); // update covariance matrix
                    

                    // calculate overall likelihood
                    likelihood = EMComputeLikelihood(this.ModelPoints, this.Segments);
                    if (Math.Abs(likelihood - prevLikelihood) < likelihoodCutoff && iterationCount > 1) { break; }
                    prevLikelihood = likelihood;
                }
            }
            catch (Exception ex)
            {
                Console.Error.WriteLine("Error fitting Gaussian mixture model: {0}", ex.Message);
                foreach (var modelPoint in this.ModelPoints)
                {
                    modelPoint.Ploidy.Mu = null;
                    modelPoint.Ploidy.Sigma = null;
                }
            }
            return likelihood;
        }
    }
}
