using System;

namespace Illumina.Common
{
    public class CrossPlatform
    {
        public static bool IsThisMono()
        {
            return Type.GetType("Mono.Runtime") != null;
        }
    }
}