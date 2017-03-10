using System.Diagnostics;
using System.Reflection;

namespace CanvasCommon
{
    public class CanvasVersionInfo
    {
        public const string NameString = "Canvas";
        public static string VersionString => Assembly.GetEntryAssembly().GetCustomAttribute<AssemblyInformationalVersionAttribute>().InformationalVersion;
        public static string CopyrightString => Assembly.GetEntryAssembly().GetCustomAttribute<AssemblyCopyrightAttribute>().Copyright;
    }
}
