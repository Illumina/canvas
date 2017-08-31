using Isas.Framework.Logging;

namespace CanvasCommon.CommandLineParsing
{
    public static class LoggerExtensions
    {
        public delegate void WriteLine(string message);
        public static WriteLine StandardWriteLine(this ILogger logger) => logger.Info;
        public static WriteLine ErrorWriteLine(this ILogger logger) => logger.Error;
    }
}