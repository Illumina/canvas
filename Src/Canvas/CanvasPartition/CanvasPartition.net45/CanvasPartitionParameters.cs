using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CanvasPartition
{
    class CanvasPartitionParameters
    {
        public int MaxInterBinDistInSegment { get; set; } = 1000000;
        public double MadFactor { get; set; } = 2.0;
        public double CBSalpha { get; set; } = 0.01;
        public double EvennessScoreThreshold { get; set; } = 94.5;
        public int EvennessScoreWindow { get; set; } = 10000;
        public double ThresholdLowerMaf { get; set; } = 0.05;
    }
}
