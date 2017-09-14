using System;
using System.Collections.Generic;
using System.Linq;
using Illumina.Common;

namespace CanvasCommon
{
    ///<summary>
    /// This class implements Density Clustering algorithm introduced in 
    /// Rodriguez, Alex, and Alessandro Laio. "Clustering by fast search and find of density peaks." Science 344.6191 (2014): 1492-1496.
    /// The principle class members are DeltaList (Delta in paper) and RhoList that are defined as follows:
    /// given distance matric d[i,j] and distanceThreshold, for each data point i compure 
    /// Densities(i) = total number of data points within distanceThreshold
    /// DistanceToHeavierNeighbor(i) = distance to nearest neighbor with higher density (min(d[i,j] for all j:=Densities(j)>Densities(i)))
    ///</summary>
    public class DensityClusteringModel
    {
        #region Members

        public List<SegmentInfo> Segments;
        private double[] _distanceArray;
        private List<double> _distanceToNearestHeavierNeighbor;
        private List<double> _densities;
        private List<double> _centroidsMafs;
        private List<double> _centroidsCoverage;
        private readonly double _coverageWeightingFactor;
        private readonly double _outlierDistanceCutoff;
        private readonly double _distanceToNearestHeavierNeighborCutoff; // used to define the centroids
        
        // parameters
        // RhoCutoff and CentroidsCutoff estimated from running density clustering on 70 HapMix tumour samples https://git.illumina.com/Bioinformatics/HapMix/
        // and visually inspecting validity of clusters
        public const double DensityCutoff = 2.0;
        public const int MaxClusterNumber = 7; // too many clusters suggest incorrect cluster partitioning 
        private const double NeighborRate = 0.02; // 1-2% is the range of values suggested in the paper

        #endregion

        private DensityClusteringModel(List<SegmentInfo> segments, double coverageWeightingFactor, double outlierDistanceCutoff, double distanceToHeavierNeighborCutoff)
        {
            Segments = segments;
            _coverageWeightingFactor = coverageWeightingFactor;
            _outlierDistanceCutoff = outlierDistanceCutoff;
            _distanceToNearestHeavierNeighborCutoff = distanceToHeavierNeighborCutoff;
        }

        public static DensityClusteringModel RunDensityClustering(List<SegmentInfo> segments, double coverageWeightingFactor, double outlierDistanceCutoff, double distanceToNearestHeavierNeighborCutoff, out int clusterCount, double densityCutoff = DensityCutoff)
        {
            // use "PloidyInfo.OutlierClusterFlag" as the default ClusterId
            segments.Select(x =>
            {
                x.ClusterId = PloidyInfo.OutlierClusterFlag;
                return x;
            });
            // Find clusters using only segments wit MAF
            var filteredDensityClustering = new DensityClusteringModel(segments.Where(segment => segment.MAF >= 0).ToList(), coverageWeightingFactor, outlierDistanceCutoff, distanceToNearestHeavierNeighborCutoff);
            filteredDensityClustering.CoverageScaling();
            filteredDensityClustering.GetDistanceArray();
            var distanceThreshold = filteredDensityClustering.EstimateDc();
            filteredDensityClustering.GaussianLocalDensity(distanceThreshold);
            var nearestNeighborHigherRho = filteredDensityClustering.CalculateDistanceToNearestHeavierNeighbor();
            clusterCount = filteredDensityClustering.FindClusters(nearestNeighborHigherRho, densityCutoff);
            // return DensityClusteringModel including segments with their ClusterId updated
            return new DensityClusteringModel(segments, coverageWeightingFactor, outlierDistanceCutoff, distanceToNearestHeavierNeighborCutoff);
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
        private void GetDistanceArray()
        {
            _distanceArray = new double[Segments.Count * (Segments.Count - 1) / 2];
            for (int i = 0; i < Segments.Count; i++)
            {
                for (int j = i + 1; j < Segments.Count; j++)
                {
                    _distanceArray[IndexMappingPointsToDistance(i, j)] = CalculateEuclideanDistance(Segments[i], Segments[j]);
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
            // similar code in original version is incorrect
            return (Segments.Count * i + j - (i + 2) * (i + 1) / 2); // all indice 0-based
        }

        /// <summary>
        /// Estimate average cluster Euclidean Distance for centroids identified by densityClustering function
        /// </summary>
        public List<double> GetCentroidsEuclideanDistance(List<double> centroidsMaFs, List<double> centroidsCoverage, int nClusters)
        {
            var centroidsEuclideanDistance = new List<double>();

            for (int clusterId = 0; clusterId < nClusters; clusterId++)
            {
                List<double> tmpDistance = new List<double>();
                foreach (SegmentInfo segment in this.Segments)
                {
                    if (segment.ClusterId.HasValue && clusterId + 1 == segment.ClusterId.Value)
                    {
                        var centroid = new SegmentInfo();
                        centroid.Coverage = centroidsCoverage[clusterId];
                        centroid.MAF = centroidsMaFs[clusterId];
                        tmpDistance.Add(CalculateEuclideanDistance(segment, centroid));
                    }
                }
                centroidsEuclideanDistance.Add(tmpDistance.Average()); 
            }
            return centroidsEuclideanDistance;
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
            return _centroidsMafs;
        }

        public List<double> GetCentroidsCoverage()
        {
            return _centroidsCoverage;
        }


        /// <summary>
        /// Return the squared euclidean distance between (coverage, maf) and (coverage2, maf2) in scaled coverage/MAF space.
        /// </summary>
        private static double CalculateEuclideanDistance(SegmentInfo segment1, SegmentInfo segment2)
        {
            return Math.Sqrt(Math.Pow(segment1.Coverage - segment2.Coverage, 2) +
                Math.Pow(segment1.MAF - segment2.MAF, 2));
        }

        /// <summary>
        /// neighborRate = average of number of elements of comb per row that are less than dc minus 1 divided by size
        /// </summary>
        private double EstimateDc(double neighborRate = NeighborRate)
        {
            if (_distanceArray.Length == 0)
                throw new Exception("Empty Distance Array!");
            // convert _distanceArray to float[] as GetPercentileNoNaNs method not implemented for double[]
            return MathSupportFunctions.GetPercentileNoNaNs(_distanceArray.Select(x => (float)x).ToArray(), 0, _distanceArray.Length - 1, (Decimal) (1 - neighborRate));
        }

        // SK: this is the method used in the paper but not used in CANVAS?
        public void NonGaussianLocalDensity(double distanceThreshold)
        {
            int ncol = Segments.Count;
            int nrow = Segments.Count;
            _densities = new List<double>(nrow);
            for (int iRho = 0; iRho < nrow; iRho++)
                this._densities.Add(0);
            int i = 0;
            for (int col = 0; col < ncol; col++)
            {
                for (int row = col + 1; row < nrow; row++)
                {
                    if (this._distanceArray[i] < distanceThreshold)
                    {
                        this._densities[row] += 1;
                        this._densities[col] += 1;
                    }
                    i++;
                }
            }
        }

        public void GaussianLocalDensity(double distanceThreshold)
        {
            int distanceLength = _distanceArray.Length;
            var half = new List<double>(new double[distanceLength]);

            for (int index = 0; index < distanceLength; index++)
            {
                double combOver = _distanceArray[index] / distanceThreshold;
                double negSq = Math.Pow(combOver, 2) * -1;
                half[index] = Math.Exp(negSq);
            }

            _densities = new List<double>(new double[Segments.Count]);
            for (int i = 0; i < Segments.Count; i++)
            {
                for (int j = i + 1; j < Segments.Count; j++)
                {
                    double temp = half[IndexMappingPointsToDistance(i, j)];
                    _densities[j] += temp;
                    _densities[i] += temp;
                }
            }
        }

        /// <summary>
        /// Estimate DistanceToNearestHeavierNeighbor value as
        /// DistanceToNearestHeavierNeighbor(i) = distance of the nearest data point of	higher density (min(d[i,j] for all j:=Densities(j)>Densities(i)))
        /// </summary>
        public Dictionary<int, int> CalculateDistanceToNearestHeavierNeighbor()
        {
            // default DistanceToNearestHeavierNeighbor value for each segment
            _distanceToNearestHeavierNeighbor = Enumerable.Repeat(double.MaxValue, Segments.Count).ToList();
            // the index of nearest segment with higher density value
            var nearestNeighborHigherRho = new Dictionary<int, int>();
            // index of segments with the highest density
            int maxDensityIndex = _densities.IndexOf(_densities.Max());
            // DistanceToNearestHeavierNeighbor[maxDensityIndex] = max distance between it and other points
            _distanceToNearestHeavierNeighbor[maxDensityIndex] = Enumerable.Range(0, Segments.Count)
                .Where(x => x != maxDensityIndex)
                .Select(x => IndexMappingPointsToDistance(x, maxDensityIndex))
                .Select(x => _distanceArray[x])
                .Max();
            nearestNeighborHigherRho.Add(maxDensityIndex, -1); // alredy the one with highest density
            for (int i = 0; i < Segments.Count; i++)
            {
                for (int j = i + 1; j < Segments.Count; j++)
                {
                    double distance = _distanceArray[IndexMappingPointsToDistance(i, j)];
                    (int indexLowerRho, int indexHigherRho) = _densities[i] < _densities[j] ? (i, j) : (j, i);
                    if (distance >= _distanceToNearestHeavierNeighbor[indexLowerRho]) continue;
                    _distanceToNearestHeavierNeighbor[indexLowerRho] = distance;
                    nearestNeighborHigherRho[indexLowerRho] = indexHigherRho;
                }
            }
            return nearestNeighborHigherRho;
        }

        private int FindClusters(Dictionary<int, int> nearestNeighborHigherRho, double rhoCutoff = DensityCutoff)
        {
            // the indice of centroids
            List<int> centroidsIndex = Enumerable.Range(0, Segments.Count)
                .Where(i => _densities[i] > rhoCutoff && _distanceToNearestHeavierNeighbor[i] > _distanceToNearestHeavierNeighborCutoff)
                .ToList();

            if (!centroidsIndex.Any())
            {
                Segments = Segments.Select(x =>
                {
                    x.ClusterId = PloidyInfo.OutlierClusterFlag;
                    return x;
                }).ToList();
                return 0;
            }

            _centroidsMafs = centroidsIndex.Select(x => Segments[x].MAF).ToList();
            _centroidsCoverage = centroidsIndex.Select(x => Segments[x].Coverage).ToList();

            // sort list and return indices
            var sortedRhoValues = _densities.Select((x, i) => new KeyValuePair<int, double>(i, x)).OrderByDescending(x => x.Value).ToList();
            var segmentIndexSortedByRho = sortedRhoValues.Select(x => x.Key).ToList();
            int clusterId = 1;
            foreach (int index in segmentIndexSortedByRho)
            {
                // set segment cluster value to the cluster centroid 
                if (centroidsIndex.Contains(index))
                {
                    Segments[index].ClusterId = clusterId;
                    clusterId++;
                    //this.Segments[index].ClusterId = CentroidsIndex.FindIndex(x => x == index) + 1;
                }
                else if (Segments[index].SumDistToKNearestNeighbours > _outlierDistanceCutoff)
                {
                    Segments[index].ClusterId = PloidyInfo.OutlierClusterFlag;
                }
                // set segment cluster value to the closest cluster segment 
                else
                {
                    // SK: could it be possible that the nearest neighbor is marked as an outlier?
                    Segments[index].ClusterId = Segments[nearestNeighborHigherRho[index]].ClusterId;
                }
            }
            return centroidsIndex.Count;
        }
    }
}