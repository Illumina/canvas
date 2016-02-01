using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Illumina.SecondaryAnalysis
{
    public interface ILogger : IDisposable
    {
        void Log(LogEntry entry);
    }

    public static class LoggerExtensions
    {
        /// <summary>
        /// Logs the specified message with Debug severity.
        /// </summary>
        /// <param name="message">The message.</param>
        public static void Debug(this ILogger logger, string message)
        {
            logger.Log(new LogEntry(LogEventType.Debug,
                message, null, null));
        }

        /// <summary>
        /// Logs the specified message with Debug severity.
        /// </summary>
        /// <param name="format">The message or format template.</param>
        /// <param name="args">Any arguments required for the format template.</param>
        public static void Debug(this ILogger logger, string format, params object[] args)
        {
            string message = string.Format(format, args);
            logger.Debug(message);
        }

        /// <summary>
        /// Logs the specified message with Info severity.
        /// </summary>
        /// <param name="message">The message.</param>
        public static void Info(this ILogger logger, string message)
        {
            logger.Log(new LogEntry(LogEventType.Information,
                message, null, null));
        }

        /// <summary>
        /// Logs the specified message with Info severity.
        /// </summary>
        /// <param name="format">The message or format template.</param>
        /// <param name="args">Any arguments required for the format template.</param>
        public static void Info(this ILogger logger, string format, params object[] args)
        {
            string message = string.Format(format, args);
            logger.Info(message);
        }

        /// <summary>
        /// Logs the specified message with Warn severity.
        /// </summary>
        /// <param name="message">The message.</param>
        public static void Warn(this ILogger logger, string message)
        {
            logger.Log(new LogEntry(LogEventType.Warning,
                message, null, null));
        }

        /// <summary>
        /// Logs the specified message with Warn severity.
        /// </summary>
        /// <param name="format">The message or format template.</param>
        /// <param name="args">Any arguments required for the format template.</param>
        public static void Warn(this ILogger logger, string format, params object[] args)
        {
            string message = string.Format(format, args);
            logger.Warn(message);
        }

        /// <summary>
        /// Logs the specified message with Error severity.
        /// </summary>
        /// <param name="message">The message.</param>
        public static void Error(this ILogger logger, string message)
        {
            logger.Log(new LogEntry(LogEventType.Error,
                message, null, null));
        }

        /// <summary>
        /// Logs the specified message with Error severity.
        /// </summary>
        /// <param name="format">The message or format template.</param>
        /// <param name="args">Any arguments required for the format template.</param>
        public static void Error(this ILogger logger, string format, params object[] args)
        {
            string message = string.Format(format, args);
            logger.Error(message);
        }

        /// <summary>
        /// Logs the specified exception with Error severity.
        /// </summary>
        /// <param name="exception">The exception.</param>
        public static void Error(this ILogger logger, Exception exception)
        {
            logger.Log(new LogEntry(LogEventType.Error,
                exception.ToString(), null, exception));
        }
    }

    public class LogEntry
    {
        public readonly LogEventType Severity;
        public readonly string Message;
        public readonly string Source;
        public readonly Exception Exception;

        /// <summary>
        /// Create a log message
        /// </summary>
        /// <param name="severity">Controls if and where the message is logged</param>
        /// <param name="message"></param>
        /// <param name="source"></param>
        /// <param name="exception"></param>
        /// <param name="addTag">Prepend severity info to message? Default if severity >= Warning</param>
        public LogEntry(LogEventType severity, string message, string source,
            Exception exception, bool? addTag=null)
        {
            if (string.IsNullOrEmpty(message))
                throw new ArgumentException($"{nameof(message)} cannot be null or empty");
            if (severity < LogEventType.Debug || severity > LogEventType.Critical)
                throw new ArgumentOutOfRangeException(nameof(severity));

            this.Severity = severity;
            if (addTag ?? (severity >= LogEventType.Warning))
                this.Message = $"{severity.ToString().ToUpper()}: {message}";
            else
                this.Message = message;
            this.Source = source;
            this.Exception = exception;
        }
    }

    public enum LogEventType
    {
        Debug = 0,      //	A debug event. This indicates a verbose event, useful during development.
        Information = 1,//	An information event. This indicates a significant, successful operation.
        Warning = 2,    //	A warning event. This indicates a problem that is not immediately significant, but that may signify conditions that could cause future problems.	
        Error = 3,      //	An error event. This indicates a significant problem the user should know about; usually a loss of functionality or data.
        Critical = 4,   //	A critical event. This indicates a fatal error or application crash.
    }
}
