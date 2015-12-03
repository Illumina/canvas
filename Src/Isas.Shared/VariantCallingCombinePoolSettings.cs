using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Isas.Shared
{
    public class VariantCallingCombinePoolSettings
    {
        public enum CombineQScoreMethod
        {
            TakeMin = 0,
            CombinePoolsAndReCalculate,
        }

        public CombineQScoreMethod HowToCombineQScore = CombineQScoreMethod.TakeMin;
        public double MaxPercentSampleSpecificError = 100; //SampleSpecificErrorRate determined by counting how often artifacts show up in one pool only.

        public VariantCallingCombinePoolSettings Clone()
        {
            return (VariantCallingCombinePoolSettings)MemberwiseClone();
        }
    }
}
