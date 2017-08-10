using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection.Metadata.Ecma335;
using Illumina.Common;
using MathNet.Numerics.LinearAlgebra.Double;

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
        public double[] DistanceArray;
        public List<double> Centroids;
        public List<double> Rho;
        private List<SegmentInfo> _filteredSegments;
        private List<double> _centroidsMafs;
        private List<double> _centroidsCoverage;
        private readonly double _coverageWeightingFactor;
        private readonly double _knearestNeighbourCutoff;
        private readonly double _centroidsCutoff;

        // parameters
        // RhoCutoff and CentroidsCutoff estimated from running density clustering on 70 HapMix tumour samples https://git.illumina.com/Bioinformatics/HapMix/
        // and visually inspecting validity of clusters
        public const double RhoCutoff = 2.0; // SK: used for outlier detection?
        // SK: MaxClusterNumber not used
        public const int MaxClusterNumber = 7; // too many clusters suggest incorrect cluster partitioning 

        private const double NeighborRateLow = 0.02;
        private const double NeighborRateHigh = 0.03;

        #endregion

        public DensityClusteringModel(List<SegmentInfo> segments, double coverageWeightingFactor, double knearestNeighbourCutoff, double centroidsCutoff)
        {
            Segments = segments;
            _coverageWeightingFactor = coverageWeightingFactor;
            _knearestNeighbourCutoff = knearestNeighbourCutoff;
            _centroidsCutoff = centroidsCutoff;
        }

        public static DensityClusteringModel RunDensityClustering(List<SegmentInfo> segments, double coverageWeightingFactor, double knearestNeighbourCutoff, double centoridCutoff, out int clusterCount, double rhoCutoff = RhoCutoff)
        {
            var densityClustering = new DensityClusteringModel(segments, coverageWeightingFactor, knearestNeighbourCutoff, centoridCutoff);
            densityClustering.SegmentsFiltering();
            densityClustering.CoverageScaling();
            densityClustering.GenDistanceArray();
            var distanceThreshold = densityClustering.EstimateDc();
            densityClustering.GaussianLocalDensity(distanceThreshold);
            densityClustering.FindCentroids();
            clusterCount = densityClustering.FindClusters(rhoCutoff);
            return densityClustering;
        }

        /// <summary>
        /// Only use segments with non-null MAF values
        /// </summary>
        private void SegmentsFiltering()
        {
            Segments = Segments.Where(segment => segment.MAF >= 0).ToList();
        }

        /// <summary>
        /// Calculate the weigthed coverage values 
        /// </summary>
        private void CoverageScaling()
        {
            Segments = Segments.Select(x =>
            {
                x.Coverage *= _coverageWeightingFactor;
                return x;
            }).ToList();
        }

        /// <summary>
        /// Compute lower triangle of the distance matrix stored by columns in a vector.  
        /// </summary>
        private void GenDistanceArray()
        {
            DistanceArray = new double[Segments.Count * (Segments.Count - 1) / 2];
            for (int i = 0; i < Segments.Count; i++)
            {
                for (int j = i + 1; j < Segments.Count; j++)
                {
                    DistanceArray[IndexMappingPointsToDistance(i, j)] = CalEuclideanDistance(Segments[i], Segments[j]);
                }
            }
        }

        // Map the indice of two different points to the index of their distance
        private int IndexMappingPointsToDistance(int i, int j)
        {
            if (i == j) throw new Exception("The two indice are the same. Please provide different ones.");
            if (i > j) // lower triangle of the distance matrix used
            {
                int tmp = i;
                i = j;
                j = tmp;
            }
            return (i * (Segments.Count - 2) + j - 1); // all indice 0-based
        }

        /// <summary>
        /// Estimate cluster variance using cluster centroids identified by densityClustering function
        /// </summary>
        public List<double> GetCentroidsVariance(List<double> centroidsMaFs, List<double> centroidsCoverage, int nClusters)
        {
            var clusterVariance = new List<double>();

            for (int clusterID = 0; clusterID < nClusters; clusterID++)
            {
                List<double> tmpDistance = new List<double>();
                foreach (SegmentInfo segment in this.Segments)
                {
                    // SK: 1-based CluserId v.s. 0-based clusterID??
                    if (segment.ClusterId.HasValue && clusterID + 1 == segment.ClusterId.Value)
                    {
                        // SK: distane between the segment and corresponding centriod
                        var centroid = new SegmentInfo();
                        centroid.Coverage = centroidsCoverage[clusterID];
                        centroid.MAF = centroidsMaFs[clusterID];
                        tmpDistance.Add(CalEuclideanDistance(segment, centroid));
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

            foreach (SegmentInfo segment in this.Segments)
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
            return _centroidsMafs;
        }

        public List<double> GetCentroidsCoverage()
        {
            return _centroidsCoverage;
        }


        /// <summary>
        /// Return the squared euclidean distance between (coverage, maf) and (coverage2, maf2) in scaled coverage/MAF space.
        /// </summary>
        private static double CalEuclideanDistance(SegmentInfo segment1, SegmentInfo segment2)
        {
            return Math.Sqrt(Math.Pow(segment1.Coverage - segment2.Coverage, 2)+ 
                Math.Pow(segment1.MAF - segment2.MAF, 2));
        }

        /// <summary>
        /// neighborRate = average of number of elements of comb per row that are less than dc minus 1 divided by size
        /// </summary>
        private double EstimateDc(double neighborRateLow = NeighborRateLow, double neighborRateHigh = NeighborRateHigh)
        {
            if (DistanceArray.Length == 0)
                throw new Exception("Empty DistanceArray!");
            var tmpLow = DistanceArray.Where(x => x > 0).Min(); // translated from old code. Why exclude zero??
            var tmpHigh = DistanceArray.Max();

            int segmentsLength = Segments.Count;
            double distanceThreshold = 0;
            var iterations = 0;
            const int maxIterations = 100000;
            while (true)
            {
                distanceThreshold = (tmpLow + tmpHigh) / 2;
                double neighborRateTmp = DistanceArray.Where(x => x < distanceThreshold).Count();
                // this part is really confusing !!!!
                if (distanceThreshold > 0) // SK: can this value <= 0?
                    neighborRateTmp = neighborRateTmp + segmentsLength; // ?

                var neighborRate = (neighborRateTmp * 2 / segmentsLength - 1) / segmentsLength;

                if (neighborRate >= neighborRateLow && neighborRate <= neighborRateHigh)
                    break;

                if (neighborRate < neighborRateLow)
                {
                    tmpLow = distanceThreshold;
                }
                else
                {
                    tmpHigh = distanceThreshold;
                }
                iterations++;
                if (iterations > maxIterations)
                    break;
            }
            return distanceThreshold;
        }

        // SK: this is the method used in the paper but not used in CANVAS?
        public void NonGaussianLocalDensity(double distanceThreshold)
        {
            int ncol = Segments.Count;
            int nrow = Segments.Count;
            Rho = new List<double>(nrow);
            for (int iRho = 0; iRho < nrow; iRho++)
                this.Rho.Add(0);
            int i = 0;
            for (int col = 0; col < ncol; col++)
            {
                for (int row = col + 1; row < nrow; row++)
                {
                    if (this.DistanceArray[i] < distanceThreshold)
                    {
                        this.Rho[row] += 1;
                        this.Rho[col] += 1;
                    }
                    i++;
                }
            }
        }

        public void GaussianLocalDensity(double distanceThreshold)
        {
            int distanceLength = DistanceArray.Length;
            var half = new List<double>(new double[distanceLength]);

            for (int index = 0; index < distanceLength; index++)
            {
                double combOver = DistanceArray[index] / distanceThreshold;
                double negSq = Math.Pow(combOver, 2) * -1;
                half[index] = Math.Exp(negSq);
            }

            Rho = new List<double>(new double[Segments.Count]);
            for (int i = 0; i < Segments.Count; i++)
            {
                for (int j = i + 1; j < Segments.Count; j++)
                {
                    double temp = half[IndexMappingPointsToDistance(i, j)];
                    Rho[j] += temp;
                    Rho[i] += temp;
                }
            }
        }

        /// <summary>
        /// Estimate Centroids value as
        /// Centroids(i) = distance	of the closes data point of	higher density (min(d[i,j] for all j:=Rho(j)>Rho(i)))
        /// </summary>
        public void FindCentroids()
        {
            Centroids = Enumerable.Repeat(double.MaxValue, Segments.Count).ToList();
            //var maximum = new List<double>(new double[Segments.Count]);
            // index of segments with the highest density
            int maxDensityIndex = Rho.IndexOf(Rho.Max());
            // Centroids[maxDensityIndex] = max distance between it and other points
            Centroids[maxDensityIndex] = Enumerable.Range(0, Segments.Count)
                .Where(x => x != maxDensityIndex)
                .Select(x => IndexMappingPointsToDistance(x, maxDensityIndex))
                .Select(x => DistanceArray[x])
                .Max();
            for (int i = 0; i < Segments.Count; i++)
            {
                for (int j = i + 1; j < Segments.Count; j++)
                {
                    double distance = DistanceArray[IndexMappingPointsToDistance(i, j)];
                    int indexLowerRho = Rho[i] < Rho[j] ? i : j;
                    if (distance < Centroids[indexLowerRho]) Centroids[indexLowerRho] = distance;
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
                return this.DistanceArray[segmentsLength * tmpIndex - (tmpIndex * (tmpIndex + 1)) / 2 + runOrderIndex - tmpIndex - 1];
            }
            else if (tmpIndex > runOrderIndex)
            {
                return this.DistanceArray[segmentsLength * runOrderIndex - (runOrderIndex * (runOrderIndex + 1)) / 2 + tmpIndex - runOrderIndex - 1];
            }
            else
            {
                return null;
            }
        }


        public int FindClusters(double rhoCutoff = RhoCutoff)
        {
            _centroidsMafs = new List<double>();
            _centroidsCoverage = new List<double>();
            int segmentsLength = this.Segments.Count;
            List<int> CentroidsIndex = new List<int>(segmentsLength);
            for (int segmentIndex = 0; segmentIndex < segmentsLength; segmentIndex++)
            {
                if (this.Rho[segmentIndex] > rhoCutoff && this.Centroids[segmentIndex] > _centroidsCutoff && this.Segments[segmentIndex].MAF >= 0)
                {
                    CentroidsIndex.Add(segmentIndex);
                    _centroidsMafs.Add(this.Segments[segmentIndex].MAF);
                    _centroidsCoverage.Add(this.Segments[segmentIndex].Coverage);
                }
            }


            // sort list and return indices
            var sortedScores = Rho.Select((x, i) => new KeyValuePair<double, int>(x, i)).OrderByDescending(x => x.Key).ToList();
            var runOrder = sortedScores.Select(x => x.Value).ToList();

            foreach (int runOrderIndex in runOrder)
            {
                // set segment cluster value to the cluster centroid 
                if (CentroidsIndex.Contains(runOrderIndex))
                {
                    this.Segments[runOrderIndex].ClusterId = CentroidsIndex.FindIndex(x => x == runOrderIndex) + 1;
                }
                // set segment cluster value to the closest cluster segment 
                else
                {
                    double? tmpDistance = null;
                    double minDistance = Double.MaxValue;
                    int minRhoElementIndex = 0;
                    for (int tmpIndex = 0; tmpIndex < segmentsLength; tmpIndex++)
                    {
                        if (Rho[tmpIndex] > Rho[runOrderIndex] && this.Segments[tmpIndex].MAF >= 0)
                        {
                            tmpDistance = GetDistance(segmentsLength, tmpIndex, runOrderIndex);

                            if (tmpDistance.HasValue)
                            {
                                if (tmpDistance < minDistance)
                                {
                                    minRhoElementIndex = tmpIndex;
                                    minDistance = (double)tmpDistance;
                                }
                            }
                        }
                    }
                    // populate clusters
                    if (this.Segments[runOrderIndex].MAF >= 0)
                        this.Segments[runOrderIndex].ClusterId = this.Segments[minRhoElementIndex].ClusterId;
                    if (!this.Segments[runOrderIndex].ClusterId.HasValue || this.Segments[runOrderIndex].MAF < 0 ||
                        this.Segments[runOrderIndex].KnearestNeighbour > this._knearestNeighbourCutoff)
                        this.Segments[runOrderIndex].ClusterId = CanvasCommon.PloidyInfo.OutlierClusterFlag;
                }
            }
            return CentroidsIndex.Count;
        }
    }
}