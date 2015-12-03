using System;

namespace Illumina.Common
{
    public static class CrossPlatform
    {
        public static bool IsThisMono()
        {
            return Type.GetType("Mono.Runtime") != null;
        }
    }
}