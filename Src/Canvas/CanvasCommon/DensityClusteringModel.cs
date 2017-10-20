using System;
using System.Collections.Generic;
using System.Linq;
using Illumina.Common;

namespace CanvasCommon
{
    ///<summary>
    /// This class implements Density Clustering algorithm introduced in 
    /// Rodriguez, Alex, and Alessandro Laio. "Clustering by fast search and find of density peaks." Science 344.6191 (2014): 1492-1496.
    /// The principle class members are Centroids (Delta in paper) and Rho that are defined as follows:
    /// given distance matric d[i,j] and distanceThreshold, for each data point i compure 
    /// Rho(i) = total number of data points within distanceThreshold
    /// Centroids(i) = distance	of the closes data point of	higher density (min(d[i,j] for all j:=Rho(j)>Rho(i)))
    ///</summary>
    public class DensityClusteringModel
    {
        #region Members

        public List<SegmentInfo> Segments;
        public List<double?> Distance;
        public List<double> Centroids;
        public List<double> Rho;
        private List<double> CentroidsMAFs;
        private List<double> CentroidsCoverage;
        private double _coverageWeightingFactor;
        private double _knearestNeighbourCutoff;
        private double CentroidsCutoff;

        // parameters
        // RhoCutoff and CentroidsCutoff estimated from running density clustering on 70 HapMix tumour samples https://git.illumina.com/Bioinformatics/HapMix/
        // and visually inspecting validity of clusters
        public const double RhoCutoff = 2.0; // SK: used for outlier detection?
        // SK: MaxClusterNumber not used
        public const int MaxClusterNumber = 7; // too many clusters suggest incorrect cluster partitioning 

        private const double NeighborRate = 0.02;

        #endregion


        public DensityClusteringModel(List<SegmentInfo> segments, double coverageWeightingFactor, double knearestNeighbourCutoff, double centroidsCutoff)
        {
            Segments = segments;
            _coverageWeightingFactor = coverageWeightingFactor;
            _knearestNeighbourCutoff = knearestNeighbourCutoff;
            CentroidsCutoff = centroidsCutoff;
        }

        /// <summary>
        /// Only use segments with non-null MAF values
        /// </summary>
        // SK: better to call this CountSegmentsWithMaf
        public int GetSegmentsForClustering(List<SegmentInfo> segments)
        {
            return segments.Count(segment => segment.MAF >= 0);
        }

        /// <summary>
        /// Estimate cluster variance using cluster centroids identified by densityClustering function
        /// </summary>
        public List<double> GetCentroidsVariance(List<double> centroidsMAFs, List<double> centroidsCoverage, int nClusters)
        {
            List<double> clusterVariance = new List<double>();

            for (int clusterId = 0; clusterId < nClusters; clusterId++)
            {
                List<double> tmpDistance = new List<double>();
                foreach (SegmentInfo segment in Segments)
                {
                    // SK: 1-based CluserId v.s. 0-based clusterID??
                    if (segment.ClusterId.HasValue && clusterId + 1 == segment.ClusterId.Value)
                    {
                        // SK: distane between the segment and corresponding centriod
                        tmpDistance.Add(GetEuclideanDistance(segment.Coverage, centroidsCoverage[clusterId], segment.MAF, centroidsMAFs[clusterId]));
                    }
                }
                clusterVariance.Add(tmpDistance.Average()); // SK: variance is the average distance?
            }
            return clusterVariance;
        }

        /// <summary>
        /// Return cluster sizes
        /// </summary>
        public List<int> GetClustersSize(int nClusters)
        {
            List<int> clustersSize = Enumerable.Repeat(0, nClusters).ToList();

            foreach (SegmentInfo segment in Segments)
            {
                if (segment.ClusterId.HasValue && segment.ClusterId.Value > 0)
                {
                    clustersSize[segment.ClusterId.Value - 1] += 1;
                }
            }
            return clustersSize;
        }

        public List<double> GetCentroidsMaf()
        {
            return CentroidsMAFs;
        }

        public List<double> GetCentroidsCoverage()
        {
            return CentroidsCoverage;
        }


        /// <summary>
        /// Return the squared euclidean distance between (coverage, maf) and (coverage2, maf2) in scaled coverage/MAF space.
        /// https://en.wikipedia.org/wiki/Euclidean_distance
        /// </summary>
        private double GetEuclideanDistance(double coverage, double coverage2, double maf, double maf2)
        {
            double diff = (coverage - coverage2) * _coverageWeightingFactor;
            double distance = diff * diff;
            diff = maf - maf2;
            distance += diff * diff;
            return Math.Sqrt(distance);
        }

        /// <summary>
        /// neighborRate = average of number of elements of comb per row that are less than dc minus 1 divided by size
        /// </summary>

        public double EstimateDc(double neighborRate = NeighborRate)
        {
            var distances = Distance.Where(x => x.HasValue).ToArray();
            if (distances.Length == 0)
                throw new Exception("Empty Distance Array!");
            // convert distances to float[] as GetPercentileNoNaNs method not implemented for double[]
            return MathSupportFunctions.GetPercentileNoNaNs(distances.Where(x => x != null).Select(x => (float)x).ToArray(), 0, distances.Length - 1, (decimal)(1 - neighborRate));
        }
        
        public void GaussianLocalDensity(double distanceThreshold)
        {
            // SK: this needs to be improved
            int distanceLength = Distance.Count;
            List<double> half = new List<double>(distanceLength);
            for (int index = 0; index < distanceLength; index++)
                half.Add(0);
            for (int index = 0; index < distanceLength; index++)
            {
                if (Distance[index].HasValue)
                {
                    double combOver = (double)Distance[index] / distanceThreshold;
                    double negSq = Math.Pow(combOver, 2) * -1;
                    half[index] = Math.Exp(negSq);
                }
            }

            int ncol = Segments.Count;
            int nrow = Segments.Count;
            Rho = new List<double>(nrow);
            for (int iRho = 0; iRho < nrow; iRho++)
                Rho.Add(0);
            int i = 0;
            for (int col = 0; col < ncol; col++)
            {
                for (int row = col + 1; row < nrow; row++)
                {
                    double temp = half[i];
                    Rho[row] += temp;
                    Rho[col] += temp;
                    i++;
                }
            }
        }

        /// <summary>
        /// Compute lower triangle of the distance matrix stored by columns in a vector. If n is the number of observations, 
        /// then for i &lt; j &lt;= n, the dissimilarity between (column) i and (row) j is retrieved from index [n*i - i*(i+1)/2 + j-i+1]. 
        /// The length of the distance vector is n*(n-1)/2.
        /// </summary>
        public void EstimateDistance()
        {
            int segmentsLength = Segments.Count;
            Distance = new List<double?>(segmentsLength * (segmentsLength - 1) / 2);
            for (int col = 0; col < segmentsLength; col++)
            {
                for (int row = col + 1; row < segmentsLength; row++)
                {
                    double? tmp = null;
                    if (Segments[col].MAF >= 0 && Segments[row].MAF >= 0)
                        tmp = GetEuclideanDistance(Segments[col].Coverage, Segments[row].Coverage, Segments[col].MAF, Segments[row].MAF);
                    Distance.Add(tmp);
                }
            }
        }


        /// <summary>
        /// Estimate Centroids value as
        /// Centroids(i) = distance	of the closes data point of	higher density (min(d[i,j] for all j:=Rho(j)>Rho(i)))
        /// </summary>
        public void FindCentroids()
        {
            // new List<double>(new double[this.Segments.Count]); 
            // Enumerable.Repeat(0D, this.Segments.Count).ToList();
            int segmentsLength = Segments.Count;
            Centroids = new List<double>(segmentsLength);
            for (int iCentroids = 0; iCentroids < segmentsLength; iCentroids++)
                Centroids.Add(0);
            List<double> maximum = new List<double>(segmentsLength);
            for (int imaximum = 0; imaximum < segmentsLength; imaximum++)
                maximum.Add(0);
            int i = 0;
            for (int col = 0; col < segmentsLength; col++)
            {
                for (int row = col + 1; row < segmentsLength; row++)
                {
                    if (!Distance[i].HasValue) // SK: better to calculate i using col and row
                    {
                        i++;
                        continue;
                    }
                    double newValue = (double)Distance[i];
                    double rhoRow = Rho[row]; // this is very confusing. better not to use row and col
                    double rhoCol = Rho[col];

                    if (rhoRow > rhoCol)
                    {
                        double centroidsCol = Centroids[col];
                        // SK: this is used to check whether this.Centroids[col] has been assigned a value after initiation
                        // ReSharper disable once CompareOfFloatsByEqualityOperator
                        if (newValue < centroidsCol || centroidsCol == 0)
                        {
                            Centroids[col] = newValue;
                        }
                    }
                    else if (newValue > maximum[col])
                    {
                        maximum[col] = newValue;
                    }

                    if (rhoCol > rhoRow)
                    {
                        double centroidsRow = Centroids[row];
                        // ReSharper disable once CompareOfFloatsByEqualityOperator
                        if (newValue < centroidsRow || centroidsRow == 0)
                        {
                            Centroids[row] = newValue;
                        }
                    }
                    else if (newValue > maximum[row])
                    {
                        maximum[row] = newValue;
                    }
                    i++;
                }
            }
            for (int j = 0; j < segmentsLength; j++)
            {
                // ReSharper disable once CompareOfFloatsByEqualityOperator
                if (Centroids[j] == 0)
                {
                    Centroids[j] = maximum[j];
                }
            }
        }


        /// <summary>
        /// Helper method for FindClusters
        /// </summary>
        private double? GetDistance(int segmentsLength, int tmpIndex, int runOrderIndex)
        {
            if (tmpIndex < runOrderIndex)
            {
                // the dissimilarity between (column) i and j is retrieved from index [n*i - i*(i+1)/2 + j-i-1].
                return Distance[segmentsLength * tmpIndex - (tmpIndex * (tmpIndex + 1)) / 2 + runOrderIndex - tmpIndex - 1];
            }
            if (tmpIndex > runOrderIndex)
            {
                return Distance[segmentsLength * runOrderIndex - (runOrderIndex * (runOrderIndex + 1)) / 2 + tmpIndex - runOrderIndex - 1];
            }
            return null;
        }


        public int FindClusters(double rhoCutoff = RhoCutoff)
        {
            CentroidsMAFs = new List<double>();
            CentroidsCoverage = new List<double>();
            int segmentsLength = Segments.Count;
            List<int> centroidIndexes = new List<int>(segmentsLength);
            for (int segmentIndex = 0; segmentIndex < segmentsLength; segmentIndex++)
            {
                if (Rho[segmentIndex] > rhoCutoff && Centroids[segmentIndex] > CentroidsCutoff && Segments[segmentIndex].MAF >= 0)
                {
                    centroidIndexes.Add(segmentIndex);
                    CentroidsMAFs.Add(Segments[segmentIndex].MAF);
                    CentroidsCoverage.Add(Segments[segmentIndex].Coverage);
                }
            }


            // sort list and return indices
            var sortedScores = Rho.Select((x, i) => new KeyValuePair<double, int>(x, i)).OrderByDescending(x => x.Key).ToList();
            var runOrder = sortedScores.Select(x => x.Value).ToList();

            foreach (int runOrderIndex in runOrder)
            {
                // set segment cluster value to the cluster centroid 
                if (centroidIndexes.Contains(runOrderIndex))
                {
                    Segments[runOrderIndex].ClusterId = centroidIndexes.FindIndex(x => x == runOrderIndex) + 1;
                }
                // set segment cluster value to the closest cluster segment 
                else
                {
                    double minDistance = double.MaxValue;
                    int minRhoElementIndex = 0;
                    for (int tmpIndex = 0; tmpIndex < segmentsLength; tmpIndex++)
                    {
                        if (Rho[tmpIndex] > Rho[runOrderIndex] && Segments[tmpIndex].MAF >= 0)
                        {
                            var tmpDistance = GetDistance(segmentsLength, tmpIndex, runOrderIndex);

                            if (tmpDistance < minDistance)
                            {
                                minRhoElementIndex = tmpIndex;
                                minDistance = (double)tmpDistance;
                            }
                        }
                    }
                    // populate clusters
                    if (Segments[runOrderIndex].MAF >= 0)
                        Segments[runOrderIndex].ClusterId = Segments[minRhoElementIndex].ClusterId;
                    if (!Segments[runOrderIndex].ClusterId.HasValue || Segments[runOrderIndex].MAF < 0 ||
                        Segments[runOrderIndex].KnearestNeighbour > _knearestNeighbourCutoff)
                        Segments[runOrderIndex].ClusterId = PloidyInfo.OutlierClusterFlag;
                }
            }
            return centroidIndexes.Count;
        }
    }
}