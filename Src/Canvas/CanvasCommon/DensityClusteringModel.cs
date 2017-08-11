using System;
using System.Collections.Generic;
using System.Linq;

namespace CanvasCommon
{
    ///<summary>
    /// This class implements Density Clustering algorithm introduced in 
    /// Rodriguez, Alex, and Alessandro Laio. "Clustering by fast search and find of density peaks." Science 344.6191 (2014): 1492-1496.
    /// The principle class members are DeltaList (Delta in paper) and RhoList that are defined as follows:
    /// given distance matric d[i,j] and distanceThreshold, for each data point i compure 
    /// RhoList(i) = total number of data points within distanceThreshold
    /// DeltaList(i) = distance	of the closes data point of	higher density (min(d[i,j] for all j:=RhoList(j)>RhoList(i)))
    ///</summary>
    public class DensityClusteringModel
    {
        #region Members

        public List<SegmentInfo> Segments;
        private double[] _distanceArray;
        private List<double> _deltaList;
        private List<double> _rhoList;
        private List<double> _centroidsMafs;
        private List<double> _centroidsCoverage;
        private readonly double _coverageWeightingFactor;
        private readonly double _knearestNeighbourCutoff;
        private readonly double _deltaCutoff; // used to define the centroids

        // parameters
        // RhoCutoff and CentroidsCutoff estimated from running density clustering on 70 HapMix tumour samples https://git.illumina.com/Bioinformatics/HapMix/
        // and visually inspecting validity of clusters
        public const double RhoCutoff = 2.0; // SK: used for outlier detection?
        // SK: MaxClusterNumber not used
        public const int MaxClusterNumber = 7; // too many clusters suggest incorrect cluster partitioning 

        private const double NeighborRateLow = 0.02;
        private const double NeighborRateHigh = 0.03;

        #endregion

        public DensityClusteringModel(List<SegmentInfo> segments, double coverageWeightingFactor, double knearestNeighbourCutoff, double deltaCutoff)
        {
            Segments = segments;
            _coverageWeightingFactor = coverageWeightingFactor;
            _knearestNeighbourCutoff = knearestNeighbourCutoff;
            _deltaCutoff = deltaCutoff;
        }

        public static DensityClusteringModel RunDensityClustering(List<SegmentInfo> segments, double coverageWeightingFactor, double knearestNeighbourCutoff, double deltaCutoff, out int clusterCount, double rhoCutoff = RhoCutoff)
        {
            var densityClustering = new DensityClusteringModel(segments, coverageWeightingFactor, knearestNeighbourCutoff, deltaCutoff);
            densityClustering.SetDefaultClusterId(PloidyInfo.OutlierClusterFlag);
            // Find clusters using only segments wit MAF
            var filteredModelWithIndex = densityClustering.GetFilteredModelWithIndex();
            var filteredDensityClustering = filteredModelWithIndex.Item1;
            var indexFilteredModel = filteredModelWithIndex.Item2;
            filteredDensityClustering.CoverageScaling();
            filteredDensityClustering.GenDistanceArray();
            var distanceThreshold = filteredDensityClustering.EstimateDc();
            filteredDensityClustering.GaussianLocalDensity(distanceThreshold);
            var nearestNeighborHigherRho = filteredDensityClustering.CalDelta();
            clusterCount = filteredDensityClustering.FindClusters(nearestNeighborHigherRho, rhoCutoff);
            // Update original model
            densityClustering.ModelUpdate(filteredDensityClustering, indexFilteredModel);
            return densityClustering;
        }

        private void SetDefaultClusterId(int? defaultClusterId)
        {
            Segments = Segments.Select(x =>
            {
                x.ClusterId = defaultClusterId;
                return x;
            }).ToList();
        }

        private void ModelUpdate(DensityClusteringModel filteredModel, List<int> indexFiltered)
        {
            foreach (int i in Enumerable.Range(0, filteredModel.Segments.Count))
            {
                int indexInOrigModel = indexFiltered[i];
                Segments[indexInOrigModel].ClusterId = filteredModel.Segments[i].ClusterId;
            }
        }

        /// <summary>
        /// Only use segments with non-null MAF values
        /// </summary>
        private Tuple<DensityClusteringModel, List<int>> GetFilteredModelWithIndex()
        {
            var filteredModel = new DensityClusteringModel(GenSegmentsWithMaf(), _coverageWeightingFactor, _knearestNeighbourCutoff, _deltaCutoff);
            return new Tuple<DensityClusteringModel, List<int>>(filteredModel, GenIndexSegmentsWithMaf());
        }

        private List<SegmentInfo> GenSegmentsWithMaf()
        {
            return Segments.Where(segment => segment.MAF >= 0).ToList();
        }

        private List<int> GenIndexSegmentsWithMaf()
        {
            return Enumerable.Range(0, Segments.Count).Where(i => Segments[i].MAF >= 0).ToList();
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
            _distanceArray = new double[Segments.Count * (Segments.Count - 1) / 2];
            for (int i = 0; i < Segments.Count; i++)
            {
                for (int j = i + 1; j < Segments.Count; j++)
                {
                    _distanceArray[IndexMappingPointsToDistance(i, j)] = CalEuclideanDistance(Segments[i], Segments[j]);
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
            return (Segments.Count * i + j - (i + 2) * (i + 1) / 2); // all indice 0-based
        }

        /// <summary>
        /// Estimate cluster variance using cluster centroids identified by densityClustering function
        /// </summary>
        public List<double> GetCentroidsVariance(List<double> centroidsMaFs, List<double> centroidsCoverage, int nClusters)
        {
            var clusterVariance = new List<double>();

            for (int clusterId = 0; clusterId < nClusters; clusterId++)
            {
                List<double> tmpDistance = new List<double>();
                foreach (SegmentInfo segment in this.Segments)
                {
                    // SK: 1-based CluserId v.s. 0-based clusterID??
                    if (segment.ClusterId.HasValue && clusterId + 1 == segment.ClusterId.Value)
                    {
                        // SK: distane between the segment and corresponding centriod
                        var centroid = new SegmentInfo();
                        centroid.Coverage = centroidsCoverage[clusterId];
                        centroid.MAF = centroidsMaFs[clusterId];
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
        private static double CalEuclideanDistance(SegmentInfo segment1, SegmentInfo segment2)
        {
            return Math.Sqrt(Math.Pow(segment1.Coverage - segment2.Coverage, 2) +
                Math.Pow(segment1.MAF - segment2.MAF, 2));
        }

        /// <summary>
        /// neighborRate = average of number of elements of comb per row that are less than dc minus 1 divided by size
        /// </summary>
        private double EstimateDc(double neighborRateLow = NeighborRateLow, double neighborRateHigh = NeighborRateHigh)
        {
            if (_distanceArray.Length == 0)
                throw new Exception("Empty DistanceArray!");
            var tmpLow = _distanceArray.Where(x => x > 0).Min(); // translated from old code. Why exclude zero??
            var tmpHigh = _distanceArray.Max();

            int segmentsLength = Segments.Count;
            double distanceThreshold = 0;
            var iterations = 0;
            const int maxIterations = 100000;
            while (true)
            {
                distanceThreshold = (tmpLow + tmpHigh) / 2;
                double neighborRateTmp = _distanceArray.Where(x => x < distanceThreshold).Count();
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
            _rhoList = new List<double>(nrow);
            for (int iRho = 0; iRho < nrow; iRho++)
                this._rhoList.Add(0);
            int i = 0;
            for (int col = 0; col < ncol; col++)
            {
                for (int row = col + 1; row < nrow; row++)
                {
                    if (this._distanceArray[i] < distanceThreshold)
                    {
                        this._rhoList[row] += 1;
                        this._rhoList[col] += 1;
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

            _rhoList = new List<double>(new double[Segments.Count]);
            for (int i = 0; i < Segments.Count; i++)
            {
                for (int j = i + 1; j < Segments.Count; j++)
                {
                    double temp = half[IndexMappingPointsToDistance(i, j)];
                    _rhoList[j] += temp;
                    _rhoList[i] += temp;
                }
            }
        }

        /// <summary>
        /// Estimate DeltaList value as
        /// DeltaList(i) = distance	of the closes data point of	higher density (min(d[i,j] for all j:=RhoList(j)>RhoList(i)))
        /// </summary>
        public Dictionary<int, int> CalDelta()
        {
            // Delta value for each segment
            _deltaList = Enumerable.Repeat(double.MaxValue, Segments.Count).ToList();
            // the index of nearest segment with higher density value
            var nearestNeighborHigherRho = new Dictionary<int, int>();
            // index of segments with the highest density
            int maxDensityIndex = _rhoList.IndexOf(_rhoList.Max());
            // DeltaList[maxDensityIndex] = max distance between it and other points
            _deltaList[maxDensityIndex] = Enumerable.Range(0, Segments.Count)
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
                    (int indexLowerRho, int indexHigherRho) = _rhoList[i] < _rhoList[j] ? (i, j) : (j, i);
                    if (distance >= _deltaList[indexLowerRho]) continue;
                    _deltaList[indexLowerRho] = distance;
                    nearestNeighborHigherRho[indexLowerRho] = indexHigherRho;
                }
            }
            return nearestNeighborHigherRho;
        }

        private int FindClusters(Dictionary<int, int> nearestNeighborHigherRho, double rhoCutoff = RhoCutoff)
        {
            // the indice of centroids
            List<int> centroidsIndex = Enumerable.Range(0, Segments.Count)
                .Where(i => _rhoList[i] > rhoCutoff && _deltaList[i] > _deltaCutoff)
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
            var sortedRhoValues = _rhoList.Select((x, i) => new KeyValuePair<int, double>(i, x)).OrderByDescending(x => x.Value).ToList();
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
                else if (Segments[index].KnearestNeighbour > _knearestNeighbourCutoff)
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