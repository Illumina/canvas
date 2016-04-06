using System.Collections.Generic;
using System.Linq;

namespace Isas.Shared
{ 
    public static class ShellExtensions
    { 
        public static char ShellQuote => Utilities.IsThisMono() ? '\'' : '"';
        public static string ShellCommandSwitch => Utilities.IsThisMono() ? "-o pipefail -c" : "/C";
        public static string ShellTab => Utilities.IsThisMono() ? "\t" : "\\t";
        public static string ShellExecutable => Utilities.IsThisMono() ? "bash" : "cmd";

        public static string WrapWithShellQuote(this IFileLocation s, string quoteString = null)
        {
            if (quoteString == null) quoteString = ShellQuote.ToString();
            return quoteString + s.FullName + quoteString;
        }

        public static string WrapWithShellQuote(this IDirectoryLocation s, string quoteString = null)
        {
            if (quoteString == null) quoteString = ShellQuote.ToString();
            return quoteString + s.FullName + quoteString;
        }

        public static string WrapWithShellQuote(this string s, string quoteString = null)
        {
            if (quoteString == null) quoteString = ShellQuote.ToString();
            return quoteString + s + quoteString;
        }

        public static string WrapWithShellQuote(this IEnumerable<string> cmd, string quoteString = null)
        {
            return string.Join(" ", cmd.Select(arg => WrapWithShellQuote(arg, quoteString)));
        }

        public static string MakeShellCommand(this string cmd)
        {
            return ShellCommandSwitch + " \" " + cmd + " \" ";
        }
    }
}
