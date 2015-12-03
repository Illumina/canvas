using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Runtime.InteropServices;

namespace Illumina.SecondaryAnalysis
{
    public static class MachineInfo
    {
        private const int ErrorInsufficientBuffer = 122;
        private const double NumBytesInGB = 1073741824.0;
// ReSharper disable InconsistentNaming
        private const double NumKBInGB = 1048576.0;
// ReSharper restore InconsistentNaming

        /// <summary>
        ///     Gets <b>true</b> if this process is running in a 64 bit
        ///     environment, <b>false</b> otherwise.
        /// </summary>
        public static bool Is64BitProcess
        {
            get { return Marshal.SizeOf(typeof(IntPtr)) == 8; }
        }

        /// <summary>
        ///     Gets <b>true</b> if this is a 64 bit Windows.
        /// </summary>
        public static bool Is64BitWindows
        {
            get
            {
                // The purpose is to know if we're running in pure 32-bit 
                // or if we're running in an emulated 32-bit environment.
                // Earlier versions of this method checked for the existence 
                // of the HKLM\Software\Wow6432Node node, but this turned 
                // out to be not realiable. Apparently, this node can exist 
                // on a 32-bit Windows as well.
                try
                {
                    string sArchitecture = Environment.GetEnvironmentVariable(
                        "PROCESSOR_ARCHITECTURE", EnvironmentVariableTarget.Machine);
                    if (sArchitecture == null)
                    {
                        return false;
                    }
                    return sArchitecture.Contains("64");
                }
                catch (NotSupportedException)
                {
                    return false;
                }
                catch (ArgumentException)
                {
                    return false;
                }
            }
        }

        /// <summary>
        ///     Returns <b>true</b> if this is a 32-bit process
        ///     running on a 64-bit server.
        /// </summary>
        public static bool IsWow64Process
        {
            get { return Is64BitWindows && !Is64BitProcess; }
        }

        [DllImport("kernel32.dll", SetLastError = true)]
        private static extern bool GetLogicalProcessorInformation(
            [Out] SystemLogicalProcessorInformatioNx86[] infos,
            ref uint infoSize);

        [return: MarshalAs(UnmanagedType.Bool)]
        [DllImport("kernel32.dll", CharSet = CharSet.Auto, SetLastError = true)]
        private static extern bool GlobalMemoryStatusEx([In, Out] MemoryStatusEx lpBuffer);

        [DllImport("kernel32.dll", SetLastError = true)]
        private static extern bool GetLogicalProcessorInformation(
            [Out] SystemLogicalProcessorInformatioNx64[] infos,
            ref uint infoSize);

        private static List<ProcessorInfo> GetProcessorInfo86()
        {
            // First we're going to execute GetLogicalProcessorInformation 
            // once to make sure that we determine the size of the data 
            // that it is going to return.
            // This call should fail with error ERROR_INSUFFICIENT_BUFFER.
            uint iReturnLength = 0;
            SystemLogicalProcessorInformatioNx86[] oDummy = null;
            bool bResult = GetLogicalProcessorInformation(oDummy,
                                                          ref iReturnLength);
            if (bResult)
            {
                throw Fail("GetLogicalProcessorInformation failed.", "x86");
            }

            // Making sure that the error code that we got back isn't that 
            // there is insufficient space in the buffer.
            int iError = Marshal.GetLastWin32Error();
            if (iError != ErrorInsufficientBuffer)
            {
                throw Fail(
                    "Insufficient space in the buffer.",
                    "x86", iError.ToString());
            }

            // Now that we know how much space we should reserve, 
            // we're going to reserve it and call 
            // GetLogicalProcessorInformation again.
            uint iBaseSize = (uint)Marshal.SizeOf(
                typeof(SystemLogicalProcessorInformatioNx86));
            uint iNumberOfElements = iReturnLength / iBaseSize;
            SystemLogicalProcessorInformatioNx86[] oData =
                new SystemLogicalProcessorInformatioNx86[iNumberOfElements];
            uint iAllocatedSize = iNumberOfElements * iBaseSize;
            if (!GetLogicalProcessorInformation(oData, ref iAllocatedSize))
            {
                throw Fail(
                    "GetLogicalProcessorInformation failed",
                    "x86",
                    Marshal.GetLastWin32Error().ToString());
            }

            // Converting the data to a list that we can easily interpret.
            List<ProcessorInfo> oList = new List<ProcessorInfo>();
            foreach (SystemLogicalProcessorInformatioNx86 oInfo in oData)
            {
                oList.Add(new ProcessorInfo(oInfo.Relationship,
                                            oInfo.Flags,
                                            oInfo.ProcessorMask));
            }
            return oList;
        }

        private static List<ProcessorInfo> GetProcessorInfo64()
        {
            // First we're going to execute GetLogicalProcessorInformation 
            // once to make sure that we determine the size of the data 
            // that it is going to return.
            // This call should fail with error ERROR_INSUFFICIENT_BUFFER.
            uint iReturnLength = 0;
            SystemLogicalProcessorInformatioNx64[] oDummy = null;
            bool bResult = GetLogicalProcessorInformation(oDummy,
                                                          ref iReturnLength);
            if (bResult)
            {
                throw Fail("GetLogicalProcessorInformation failed.", "x64");
            }

            // Making sure that the error code that we got back is not  
            // that there is in sufficient space in the buffer.
            int iError = Marshal.GetLastWin32Error();
            if (iError != ErrorInsufficientBuffer)
            {
                throw Fail(
                    "Insufficient space in the buffer.",
                    "x64", iError.ToString());
            }

            // Now that we know how much space we should reserve, 
            // we're going to reserve it and call 
            // GetLogicalProcessorInformation again.
            uint iBaseSize = (uint)Marshal.SizeOf(
                typeof(SystemLogicalProcessorInformatioNx64));
            uint iNumberOfElements = iReturnLength / iBaseSize;
            SystemLogicalProcessorInformatioNx64[] oData =
                new SystemLogicalProcessorInformatioNx64[iNumberOfElements];
            uint iAllocatedSize = iNumberOfElements * iBaseSize;
            if (!GetLogicalProcessorInformation(oData, ref iAllocatedSize))
            {
                throw Fail("GetLogicalProcessorInformation failed",
                           "x64", Marshal.GetLastWin32Error().ToString());
            }

            // Converting the data to a list that we can easily interpret.
            List<ProcessorInfo> oList = new List<ProcessorInfo>();
            foreach (SystemLogicalProcessorInformatioNx64 oInfo in oData)
            {
                oList.Add(new ProcessorInfo(
                              oInfo.Relationship,
                              oInfo.Flags,
                              oInfo.ProcessorMask));
            }
            return oList;
        }

        private static Exception Fail(params string[] data)
        {
            return new NotSupportedException(
                "GetPhysicalProcessorCount unexpectedly failed " +
                "(" + String.Join(", ", data) + ")");
        }

        /// <summary>
        ///     returns the number of GB installed on the machine
        /// </summary>
        public static double TotalPhysicalMemoryGB()
        {
            if (Utilities.IsThisMono())
            {
                PerformanceCounter counter = new PerformanceCounter("Mono Memory", "Total Physical Memory");
                double physicalGB = (counter.RawValue / NumBytesInGB);
                return physicalGB;
            }
            MemoryStatusEx stat = new MemoryStatusEx();
            GlobalMemoryStatusEx(stat);
            double result = stat.TotalPhysical / NumBytesInGB;
            return result;
        }

        /// <summary>
        ///     returns the number of GB available on the machine
        /// </summary>
        public static double TotalPhysicalAvailableMemoryGB()
        {
            if (Utilities.IsThisMono())
            {
                UInt64 freeMemoryKB = 0;
                UInt64 cachedMemoryKB = 0;
                UInt64 bufferedMemoryKB = 0;
                try
                {
                    // If we're on linux, go looking for memory info
                    // in the meminfo file.
                    // Hoping to see an entry that looks like
                    //  MemFree:       1976408 kB
                    //  Buffers:        348648 kB
                    //  Cached:        4665328 kB
                    // that specifies the available memory. There is certain
                    // memory taken by the operating system (buffers, cache) that
                    // can be made available if requested.
                    StreamReader reader = new StreamReader("/proc/meminfo");
                    string input = null;
                    while ((input = reader.ReadLine()) != null)
                    {
                        string[] splat = input.Split(new[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                        if (splat.Length >= 2)
                        {
                            if (splat[0] == "MemFree:")
                            {
                                UInt64.TryParse(splat[1], out freeMemoryKB);
                            }
                            if (splat[0] == "Buffers:")
                            {
                                UInt64.TryParse(splat[1], out bufferedMemoryKB);
                            }
                            if (splat[0] == "Cached:")
                            {
                                UInt64.TryParse(splat[1], out cachedMemoryKB);
                            }
                        }
                    }
                    reader.Close();
                }
                catch (Exception e)
                {
                    Console.WriteLine(e);
                }
                return (freeMemoryKB + cachedMemoryKB + bufferedMemoryKB) / NumKBInGB;
            }
            MemoryStatusEx stat = new MemoryStatusEx();
            GlobalMemoryStatusEx(stat);
            return stat.AvailablePhysical / NumBytesInGB;
        }

        /// <summary>
        ///     Returns the number of physical processor cores (this is less than the number of virtual cores if there's hyperthreading)
        /// </summary>
        public static int GetPhysicalProcessorCoreCount()
        {
            int result = 0;
            try
            {
                result = Utilities.IsThisMono() ? GetPhysicalProcessorCountMono() : GetPhysicalProcessorCountWindows();
            }
            catch
            {
                // Silently return a reasonable default.
                result = Environment.ProcessorCount;
            }
            if (result < 1) result = Environment.ProcessorCount; // still more sanity-checking
            return result;
        }

        /// <summary>
        /// Getting the physical processor count is tricky!  We'd tried getting the info from /proc/cpuinfo, but with limited success.
        /// In theory each (Phyiscal ID, Core ID) combination in /proc/cpuinfo represents a distinct physical cores; there will 
        /// be 2 occurrences of each if hyperthreading is enabled.  See:
        /// http://www.linuxforums.org/articles/finding-server-is-multi-processor-multi-core-or-hyperthreading-is-enabled-or-not-_856.html
        /// However, this logic is not always correct!  There are amazon nodes with 8 cores where each one has
        /// physical ID = 0 and core ID = 0.  For now, we'll just treat all cores equally regardless of whether hyperthreading
        /// is enabled.  See old changesets for the (non-working) logic to use /proc/cpuinfo
        /// </summary>
        private static int GetPhysicalProcessorCountMono()
        {
            return Environment.ProcessorCount;
        }

        private static int GetPhysicalProcessorCountWindows()
        {
            int processorCount = Environment.ProcessorCount;
            if (!Is64BitProcess)
            {
                Version oVersion = Environment.OSVersion.Version;
                if (oVersion < new Version(5, 1, 2600))
                {
                    return processorCount;
                }
                if (oVersion.Major == 5 && oVersion.Minor == 1 &&
                    !Environment.OSVersion.ServicePack.Equals("Service Pack 3", StringComparison.OrdinalIgnoreCase))
                {
                    return processorCount;
                }
            }

            // Getting a list of processor information
            List<ProcessorInfo> oList = Is64BitProcess ? GetProcessorInfo64() : GetProcessorInfo86();

            // The list will basically contain something like this at this point:
            //
            // E.g. for a 2 x single core
            // Relationship              Flags      ProcessorMask
            // ---------------------------------------------------------
            // RelationProcessorCore     0          1
            // RelationProcessorCore     0          2
            // RelationNumaNode          0          3
            //
            // E.g. for a 2 x dual core
            // Relationship              Flags      ProcessorMask
            // ---------------------------------------------------------
            // RelationProcessorCore     1          5
            // RelationProcessorCore     1          10
            // RelationNumaNode          0          15
            //
            // E.g. for a 1 x quad core
            // Relationship              Flags      ProcessorMask
            // ---------------------------------------------------------
            // RelationProcessorCore     1          15
            // RelationNumaNode          0          15
            //
            // E.g. for a 1 x dual core
            // Relationship              Flags      ProcessorMask  
            // ---------------------------------------------------------
            // RelationProcessorCore     0          1              
            // RelationCache             1          1              
            // RelationCache             1          1              
            // RelationProcessorPackage  0          3              
            // RelationProcessorCore     0          2              
            // RelationCache             1          2              
            // RelationCache             1          2              
            // RelationCache             2          3              
            // RelationNumaNode          0          3
            // 
            // Vista or higher will return one RelationProcessorPackage 
            // line per socket. On other operating systems we need to 
            // interpret the RelationProcessorCore lines.
            //
            // More information:
            // http://msdn2.microsoft.com/en-us/library/ms683194(VS.85).aspx
            // http://msdn2.microsoft.com/en-us/library/ms686694(VS.85).aspx

            // return the number of non-hyperthreaded cores
            int iCount = 0;
            foreach (ProcessorInfo oItem in oList)
            {
                if (oItem.Relationship ==
                    RelationProcessorCore.RelationProcessorCore)
                {
                    iCount++;
                }
            }

            if (iCount > 0)
            {
                return iCount;
            }

            throw Fail("No cpus have been detected.");
        }

        [StructLayout(LayoutKind.Sequential)]
        private struct CacheDescriptor
        {
            private readonly byte Level;
            private readonly byte Associativity;
            private readonly UInt16 LineSize;
            private readonly UInt32 Size;
            [MarshalAs(UnmanagedType.U4)]
            private readonly ProcessorCacheType Type;
        }

        [StructLayout(LayoutKind.Sequential, CharSet = CharSet.Auto)]
        private class MemoryStatusEx
        {
            private uint Length;
            public uint MemoryLoad;
            public ulong TotalPhysical;
            public ulong AvailablePhysical;
            public ulong TotalPageFile;
            public ulong AvailablePageFile;
            public ulong TotalVirtual;
            public ulong AvailableVirtual;
            public ulong AvailableExtendedVirtual;

            public MemoryStatusEx()
            {
                Length = (uint)Marshal.SizeOf(typeof(MemoryStatusEx));
            }
        }

        private enum ProcessorCacheType
        {
            /// <summary>
            ///     The cache is unified.
            /// </summary>
            UnifiedCache = 0,

            /// <summary>
            ///     InstructionThe cache is for processor instructions.
            /// </summary>
            InstructionCache = 1,

            /// <summary>
            ///     The cache is for data.
            /// </summary>
            DataCache = 2,

            /// <summary>
            ///     TraceThe cache is for traces.
            /// </summary>
            TraceCache = 3
        }

        private class ProcessorInfo
        {
            private readonly byte _flags;
            private readonly uint _processorMask;
            private readonly RelationProcessorCore _relationship;

            public ProcessorInfo(RelationProcessorCore relationShip,
                                 byte flags, uint processorMask)
            {
                _relationship = relationShip;
                _flags = flags;
                _processorMask = processorMask;
            }

            public RelationProcessorCore Relationship
            {
                get { return _relationship; }
            }

            public byte Flags
            {
                get { return _flags; }
            }

            public uint ProcessorMask
            {
                get { return _processorMask; }
            }
        }

        private enum RelationProcessorCore
        {
            /// <summary>
            ///     The specified logical processors share a
            ///     single processor core.
            /// </summary>
            RelationProcessorCore = 0,

            /// <summary>
            ///     The specified logical processors are part
            ///     of the same NUMA node.
            /// </summary>
            RelationNumaNode = 1,

            /// <summary>
            ///     The specified logical processors  share a cache.
            ///     Windows Server 2003:  This value is not supported
            ///     until Windows Server 2003 SP1 and Windows XP
            ///     Professional x64 Edition.
            /// </summary>
            RelationCache = 2,

            /// <summary>
            ///     The specified logical processors share a physical
            ///     package (a single package socketed or soldered
            ///     onto a motherboard may contain multiple processor
            ///     cores or threads, each of which is treated as a
            ///     separate processor by the operating system).
            ///     Windows Server 2003:  This value is not
            ///     supported until Windows Vista.
            /// </summary>
            RelationProcessorPackage = 3
        }

        [StructLayout(LayoutKind.Explicit)]
        private struct SystemLogicalProcessorInformatioNx64
        {
            [FieldOffset(0)]
            public readonly uint ProcessorMask;
            [FieldOffset(8), MarshalAs(UnmanagedType.U4)]
            public readonly RelationProcessorCore Relationship;
            [FieldOffset(12)]
            public readonly byte Flags;
            [FieldOffset(12)]
            private readonly CacheDescriptor Cache;
            [FieldOffset(12)]
            private readonly UInt32 NodeNumber;
            [FieldOffset(12)]
            private readonly UInt64 Reserved1;
            [FieldOffset(20)]
            private readonly UInt64 Reserved2;
        }

        [StructLayout(LayoutKind.Explicit)]
        private struct SystemLogicalProcessorInformatioNx86
        {
            [FieldOffset(0)]
            public readonly uint ProcessorMask;
            [FieldOffset(4), MarshalAs(UnmanagedType.U4)]
            public readonly RelationProcessorCore Relationship;
            [FieldOffset(8)]
            public readonly byte Flags;
            [FieldOffset(8)]
            private readonly CacheDescriptor Cache;
            [FieldOffset(8)]
            private readonly UInt32 NodeNumber;
            [FieldOffset(8)]
            private readonly UInt64 Reserved1;
            [FieldOffset(16)]
            private readonly UInt64 Reserved2;
        }
    }
}