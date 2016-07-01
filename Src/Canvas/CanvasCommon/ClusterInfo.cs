using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CanvasCommon
{

    /// <summary>
    /// Stores information about clustering of CanvasSegment objects
    /// </summary>
    public class ClusterInfo

    {
        #region Members
        public const int OutlierClusterFlag = -1;
        public const int UndersegmentedClusterFlag = -2;
        public int ClusterId;
        public bool IsHeterogeneous = false;
        public double ClusterMedianDistance;
        public double ClusterEntropy;
        public double ClusterMeanDistance;
        public double ClusterVariance;
        public List<double> ClusterDistance = new List<double>();
        public List<double> ClusterMajorChromosomeCount = new List<double>();
        #endregion

        /// <summary>
        /// Cluster entropy is estimated Based on the following logic: 
        /// for each segment in the cluster, identify the closest model point (copy number and MCC). 
        /// Clusters in which all segments belong to only one model point will have very low entropy. 
        /// Conversely clusters containing segment that are closest to different model points will have high entropy.
        /// Clusters with high entropy tend to be located between different model points and are likely to represent 
        /// sunclonal CMV variants
        /// </summary>
                
        public void ComputeClusterEntropy()
        {
                var uniqueMccCounts = ClusterMajorChromosomeCount.Distinct().ToList();
                List<int> uniqueMccCountsIndices = Enumerable.Repeat(0, uniqueMccCounts.Count).ToList();
                foreach (double mcc in ClusterMajorChromosomeCount)
                    uniqueMccCountsIndices[uniqueMccCounts.FindIndex(x => x == mcc)] += 1;
                foreach (double mccCount in uniqueMccCounts)
                {
                    if (mccCount > 0)
                    {
                        double mccTmpProbability = mccCount / ClusterMajorChromosomeCount.Count;
                    this.ClusterEntropy += -mccTmpProbability * Math.Log(mccTmpProbability);
                    }
                }
        }
    }
}
