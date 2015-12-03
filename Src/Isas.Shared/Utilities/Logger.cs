using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Threading;
using System.Threading.Tasks;
using Illumina.SecondaryAnalysis;

namespace Isas.Shared.Utilities
{
    public class Logger : ILogger
    {
        private readonly List<TextWriter> _logWriters = new List<TextWriter>();
        private readonly List<TextWriter> _errorWriters = new List<TextWriter>();
        private const int QueueLengthThreshold = 0; // clear queues if we get more than this many log entries
        private readonly StreamWriter _workflowLog;
        private readonly StreamWriter _workflowError;
        private readonly FileStream _workflowLogFileStream;
        private readonly FileStream _workflowErrorFileStream;
        private readonly LogEventType _lowestEventTypeForWorkflowLog;
        private readonly LogEventType _lowestEventTypeForWorkflowError;
        private bool _disposed;
        private readonly object _disposeLock = new object();

        private readonly ConcurrentQueue<string> _logEntries = new ConcurrentQueue<string>();
        private readonly ConcurrentQueue<string> _errorEntries = new ConcurrentQueue<string>();
        private readonly CancellationTokenSource _cancellationTokenSource = new CancellationTokenSource();
        private TaskCompletionSource<bool> _readyToFlushLog = new TaskCompletionSource<bool>();
        private TaskCompletionSource<bool> _readyToFlushError = new TaskCompletionSource<bool>();
        private readonly Task _logTask;
        private readonly Task _errorTask;

        public Logger(IFileLocation workflowLog, IFileLocation workflowError, bool alsoWriteToStandardOutput = true, LogEventType lowestEventTypeForWorkflowLog = LogEventType.Information,
            LogEventType lowestEventTypeForWorkflowError = LogEventType.Warning) : this(alsoWriteToStandardOutput, lowestEventTypeForWorkflowLog, lowestEventTypeForWorkflowError)
        {
            _logWriters.Add(InitLog(workflowLog, out _workflowLogFileStream, out _workflowLog));
            _errorWriters.Add(InitLog(workflowError, out _workflowErrorFileStream, out _workflowError));
        }

        public Logger(IEnumerable<TextWriter> logWriters, IEnumerable<TextWriter> errorWriters,
            LogEventType lowestEventTypeForWorkflowLog = LogEventType.Information,
            LogEventType lowestEventTypeForWorkflowError = LogEventType.Warning) : this(false, lowestEventTypeForWorkflowLog, lowestEventTypeForWorkflowError)
        {
            _logWriters.AddRange(logWriters);
            _errorWriters.AddRange(errorWriters);
        }

        public Logger(bool includeStandardOutput = true, LogEventType lowestEventTypeForWorkflowLog = LogEventType.Information, LogEventType lowestEventTypeForWorkflowError = LogEventType.Warning)
        {
            if (includeStandardOutput)
            {
                _logWriters.Add(Console.Out);
                _errorWriters.Add(Console.Error);
            }

            if (lowestEventTypeForWorkflowError <= lowestEventTypeForWorkflowLog)
                throw new ArgumentException("Lowest severity level for error log must be greater than lowest severity level for standard log");
            _lowestEventTypeForWorkflowLog = lowestEventTypeForWorkflowLog;
            _lowestEventTypeForWorkflowError = lowestEventTypeForWorkflowError;

            //make sure 
            _cancellationTokenSource.Token.Register(() => _readyToFlushLog.TrySetResult(true));
            _cancellationTokenSource.Token.Register(() => _readyToFlushError.TrySetResult(true));

            //start our logging tasks
            _logTask = RunLogger(true);
            _errorTask = RunLogger(false);
        }

        private static TextWriter InitLog(IFileLocation log, out FileStream fileStream, out StreamWriter streamWriter)
        {
            fileStream = new FileStream(log.FullName, FileMode.Append, FileAccess.Write, FileShare.ReadWrite);
            streamWriter = new StreamWriter(fileStream);
            return streamWriter;
        }

        private async Task RunLogger(bool isWorkflowLog)
        {
            ConcurrentQueue<string> entries = isWorkflowLog ? _logEntries : _errorEntries;
            Task readyToFlush = isWorkflowLog ? _readyToFlushLog.Task : _readyToFlushError.Task;
            List<TextWriter> writers = isWorkflowLog ? _logWriters : _errorWriters;
            while (!_cancellationTokenSource.Token.IsCancellationRequested)
            {
                await readyToFlush.ConfigureAwait(false);
                await FlushLog(writers, entries).ConfigureAwait(false);

                // create a new TaskCompletionSource for the next batch of entries
                if (isWorkflowLog)
                {
                    _readyToFlushLog = new TaskCompletionSource<bool>();
                    readyToFlush = _readyToFlushLog.Task;
                }
                else
                {
                    _readyToFlushError = new TaskCompletionSource<bool>();
                    readyToFlush = _readyToFlushError.Task;
                }
            }
            await FlushLog(writers, entries).ConfigureAwait(false);
        }

        private static async Task FlushLog(List<TextWriter> writers, ConcurrentQueue<string> entries)
        {
            string line;
            while (entries.TryDequeue(out line))
            {
                foreach (var writer in writers)
                {
                    await writer.WriteLineAsync(line).ConfigureAwait(false);

                    await writer.FlushAsync().ConfigureAwait(false);
                }
            }
        }

        public void Log(LogEntry entry)
        {
            lock (_disposeLock)
            {
                if (_disposed) throw new ObjectDisposedException(typeof(Logger).FullName);

                ConcurrentQueue<string> entries = null;
                TaskCompletionSource<bool> readyToFlush = null;
                if (entry.Severity >= _lowestEventTypeForWorkflowError)
                {
                    entries = _errorEntries;
                    readyToFlush = _readyToFlushError;
                }
                else if (entry.Severity >= _lowestEventTypeForWorkflowLog)
                {
                    entries = _logEntries;
                    readyToFlush = _readyToFlushLog;
                }

                if (entries == null) return;

                entries.Enqueue(FormatLogMessage(entry.Message));
                if (entries.Count > QueueLengthThreshold)
                    readyToFlush.TrySetResult(true);
            }
        }

        private string FormatLogMessage(string message)
        {
            return string.Format("{0},{1},{2}", DateTime.Now.ToShortDateString(), DateTime.Now.ToString("HH:mm:ss.fff"), message);
        }

        public void Dispose()
        {
            lock (_disposeLock)
            {
                if (_disposed) return;
                _cancellationTokenSource.Cancel();
                Task.WaitAll(_logTask, _errorTask);
                _cancellationTokenSource.Dispose();
                _workflowLog?.Dispose();
                _workflowError?.Dispose();
                _workflowLogFileStream?.Dispose();
                _workflowErrorFileStream?.Dispose();
                _disposed = true;
            }
        }
    }
}
