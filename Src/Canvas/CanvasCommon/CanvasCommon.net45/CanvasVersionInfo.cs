using System.Diagnostics;
using System.Reflection;

namespace CanvasCommon
{
    public class CanvasVersionInfo
    {
        public const string NameString = "Canvas";
        public static string VersionString => FileVersionInfo.GetVersionInfo(Assembly.GetExecutingAssembly().Location).ProductVersion;
        public static string CopyrightString => FileVersionInfo.GetVersionInfo(Assembly.GetExecutingAssembly().Location).LegalCopyright;
    }
}
