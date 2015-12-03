using System;
using System.Runtime.InteropServices;
using System.Text;

namespace Illumina.Win32
{
    public static class Functions
    {

        [DllImport("kernel32.dll")]
        public static extern uint GetLastError();

        [DllImport("kernel32.dll", CharSet = CharSet.Auto)]
        public static extern bool WriteFile(
            IntPtr hFile,
            IntPtr lpBuffer,
            UInt32 nNumberOfBytesToWrite,
            out UInt32 lpNumberOfBytesWritten,
            IntPtr lpOverlapped
            );

        [DllImport("kernel32.dll", CharSet = CharSet.Auto)]
        public static extern bool ReadFile(
            IntPtr hFile,
            IntPtr lpBuffer,
            UInt32 nNumberOfBytesToRead,
            out UInt32 lpNumberOfBytesRead,
            IntPtr lpOverlapped
            );

        [DllImport("kernel32.dll", CharSet = CharSet.Auto)]
        public static extern UInt32 SetFilePointer(
            IntPtr hFile,
            Int32 lDistanceToMove,
            IntPtr lpDistanceToMoveHigh, // UInt32 *, or null
            UInt32 dwMoveMethod
            );

    }
}