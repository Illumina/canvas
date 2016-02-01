using System;
using System.ComponentModel;
using System.Runtime.InteropServices;

namespace Illumina.Common
{
    public static class CrossPlatform
    {
        public static bool IsThisMono()
        {
            return Type.GetType("Mono.Runtime") != null;
        }

        public static Win32Exception GetLastWin32Exception()
        {
            return new Win32Exception(Marshal.GetLastWin32Error());
        }
    }
}