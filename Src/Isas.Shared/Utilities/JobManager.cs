using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Runtime.CompilerServices;
using System.Runtime.Serialization;
using System.Threading;

namespace Illumina.SecondaryAnalysis
{
    /// <summary>
    /// JobManager is responsible for launching external processes. It handles output redirection,
    /// process monitoring (duration and memory usage).
    /// </summary>
    public sealed class JobManager
    {
        #region members
        private ErrorHandler _error;
        private ErrorHandler _log;

        private readonly long _maximumBytesPerProcess;
        private readonly float _maximumHoursPerProcess;

        private readonly string _alternativeWorkingDirectory;

        private const int MonitorFrequencyMilliseconds = 1000;

        private Semaphore _jobPool;
        private int _jobsRemaining;
        private AutoResetEvent _doneEvent;

        private static Object _commandLogLock = new Object();
        private static Object _rUsageLock = new Object();
        #endregion

        // constructor
        public JobManager(float maximumGigabytesPerProcess, float maximumHoursPerProcess, ErrorHandler log, ErrorHandler error)
        {
            //todo Change JobManager to track total memory used by all processes, not just maxPerProcess
            _error = error;
            _log = log;
            _maximumBytesPerProcess = (long)(maximumGigabytesPerProcess * 1073741824);
            _maximumHoursPerProcess = maximumHoursPerProcess;

            //we don't want core dumps in our Isas installation directory
            _alternativeWorkingDirectory = Directory.GetCurrentDirectory();
        }

        #region process utility methods
        /// <summary>
        /// deletes the specified log file
        /// </summary>
        private void DeleteLog(string logPath)
        {
            try
            {
                if (File.Exists(logPath)) File.Delete(logPath);
            }
            catch (Exception e)
            {
                OnError(string.Format("Unable to delete the log file ({0}): {1}", logPath, e.Message));
            }
        }

        /// <summary>
        /// returns the total walltime for this process
        /// </summary>
        private static TimeSpan? GetProcessWallTime(Process jobProcess)
        {
            try
            {
                return DateTime.Now.Subtract(jobProcess.StartTime);
            }
            catch
            {
                return null; // the process has already exited
            }
        }

        /// <summary>
        /// returns the exit code for the specified process. Handles strange situations where the process is not
        /// quite finished.
        /// </summary>
        private static int GetProcessExitCode(Process jobProcess)
        {
            bool needExitCode = true;
            int exitCode = 0;

            while (needExitCode)
            {
                try
                {
                    exitCode = jobProcess.ExitCode;
                    needExitCode = false;
                }
                catch (InvalidOperationException)
                {
                    // even though we used jobProcess.WaitForExit, we sometimes get an InvalidOperationException
                    // in mono because the process has not exited yet *shrug*. So let's wait and try again.
                    Thread.Sleep(1000);
                }
            }

            return exitCode;
        }

        /// <summary>
        /// returns the number of bytes used by this process
        /// </summary>
        private static long? GetProcessMemoryUsage(Process jobProcess)
        {
            try
            {
                return jobProcess.WorkingSet64;
            }
            catch
            {
                return null; // the process has already exited
            }
        }

        /// <summary>
        /// returns the number of bytes in virtual address space used by this process
        /// </summary>
        private static long? GetPeakVirtualMemorySize64(Process jobProcess)
        {
            try
            {
                return jobProcess.PeakVirtualMemorySize64;
            }
            catch
            {
                return null; // the process has already exited
            }
        }

        /// <summary>
        /// returns true if we successfully killed the process. Otherwise the process probably already ended.
        /// </summary>
        private static bool KillProcess(Process jobProcess)
        {
            bool successfulKill = true;

            try
            {
                jobProcess.Kill();
            }
            catch (SystemException)
            {
                // this is sometimes thrown if the process has already finished.
                successfulKill = false;
            }

            return successfulKill;
        }
        #endregion


        /// <summary>
        /// Get up to 5 lines (or 2k) of logging, whichever comes first.
        /// </summary>
        private static string GetLastLines(string logPath)
        {
            if (!File.Exists(logPath)) return null;
            byte[] Buffer = new byte[1024 * 2];
            try
            {
                using (FileStream logStream = new FileStream(logPath, FileMode.Open))
                {
                    long seekPos = Math.Max(0, logStream.Length - Buffer.Length);
                    logStream.Seek(seekPos, SeekOrigin.Begin);
                    int bytesRead = logStream.Read(Buffer, 0, Buffer.Length);
                    int newlineCount = 0;
                    int newlinePos = -1;
                    for (int pos = bytesRead - 1; pos >= 0; pos--)
                    {
                        if (Buffer[pos] == '\n')
                        {
                            newlineCount++;
                            newlinePos = pos;
                            if (newlineCount >= 5) break;
                        }
                    }
                    if (newlinePos == -1) return null;
                    System.Text.Encoding encoding = new System.Text.ASCIIEncoding();
                    return encoding.GetString(Buffer, newlinePos, bytesRead - newlinePos).Trim();
                }
            }
            catch
            {
                return null;
            }
        }

        /// <summary>
        /// outputs an error message or throws an exception if a job returns a non-zero exit code
        /// </summary>
        private void CheckExitCodes(List<UnitOfWork> jobs, bool throwExceptionOnBadExitCode)
        {
            foreach (UnitOfWork job in jobs)
            {
                if (job.ExitCode != 0)
                {
                    string exceptionMessage = string.Empty;
                    if (job.ThrownException != null)
                    {
                        exceptionMessage = string.Format("{0} returned error code {1}\n** Failed command:\n{2} {3}\n** Stderr output (if any) is at: {4}\n  Full error message: {5}",
                            Path.GetFileName(job.ExecutablePath), job.ExitCode, job.ExecutablePath, job.CommandLine, job.ErrorLogPath,
                            job.ThrownException);
                    }
                    else
                    {
                        // Look for output in the stderr log file:
                        string stderr = GetLastLines(job.ErrorLogPath);

                        exceptionMessage = string.Format("{0} returned error code {1}\n\n---Failed command:\n{2} {3}",
                            Path.GetFileName(job.ExecutablePath), job.ExitCode,
                            job.ExecutablePath, job.CommandLine);
                        if (!string.IsNullOrEmpty(stderr))
                        {
                            exceptionMessage += string.Format("\n\n---Last lines of stderr log file from {0}:\n{1}\n\n", job.ErrorLogPath, stderr);
                        }
                        else
                        {
                            string stdout = GetLastLines(job.OutputLogPath);
                            if (!string.IsNullOrEmpty(stdout))
                                exceptionMessage += string.Format("\n\n--- stderr log file is empty! Last lines of stdout log file from {0}:\n{1}\n\n", job.OutputLogPath, stdout);
                            else
                                exceptionMessage += string.Format("\n\n--- both stdout and stderr log files are empty!");
                        }
                    }
                    // throw exceptions
                    if (throwExceptionOnBadExitCode)
                    {
                        throw new JobFailedException(exceptionMessage);
                    }

                    // display error message
                    OnError(exceptionMessage);
                }
            }
        }

        /// <summary>
        /// returns a new command line string where the Isas tokens have been replaced with the appropriate values
        /// </summary>
        private static string ConvertCommandLineTokens(string commandLine)
        {
            return commandLine.Replace(CommandLineTokens.NumProcessors, Environment.ProcessorCount.ToString(CultureInfo.InvariantCulture));
        }

        /// <summary>
        /// Monitors a running process to ensure that we don't use too much memory or run too long
        /// </summary>
        private void MonitorProcess(Process jobProcess, ref bool stopMonitorThread, ref long maxVMem)
        {
            try
            {
                while (true)
                {
                    // quit if the process has already finished
                    if (stopMonitorThread) break;

                    // refresh the process information
                    jobProcess.Refresh();
                    TimeSpan? wallTime = GetProcessWallTime(jobProcess);
                    long? vMem = GetProcessMemoryUsage(jobProcess);

                    // Update max RAM used by any process in this batch
                    if (GetPeakVirtualMemorySize64(jobProcess).GetValueOrDefault(0) > maxVMem)
                        maxVMem = jobProcess.PeakVirtualMemorySize64;

                    // check the wall time
                    if (wallTime.HasValue && wallTime.Value.TotalHours > _maximumHoursPerProcess)
                    {
                        if (KillProcess(jobProcess))
                        {
                            OnError(string.Format("Job Manager detected that {0} ran for more than {1:0.00} hours. Killing process.",
                                jobProcess.StartInfo.FileName,
                                _maximumHoursPerProcess));
                            break;
                        }
                    }

                    // check the memory usage
                    if (vMem.GetValueOrDefault(0) > _maximumBytesPerProcess)
                    {
                        if (KillProcess(jobProcess))
                        {
                            OnError(string.Format("Job Manager detected that {0} used more than {1:0.00} GB of memory. Killing process.",
                                jobProcess.StartInfo.FileName,
                                _maximumBytesPerProcess / (float)(1024 * 1024 * 1024)));
                            break;
                        }
                    }

                    Thread.Sleep(MonitorFrequencyMilliseconds);
                }
            }
            catch (ThreadAbortException)
            {
                if (!stopMonitorThread) KillProcess(jobProcess);
            }
        }

        /// <summary>
        ///     writes to the standard error log
        /// </summary>
        private void OnError(string message)
        {
            if (_error != null) _error(message);
            else Console.WriteLine(message);
        }

        /// <summary>
        ///     writes to the standard output log
        /// </summary>
        private void OnLog(string message)
        {
            if (_log != null) _log(message);
            else Console.WriteLine(message);
        }

        /// <summary>
        /// processes a list of jobs
        /// </summary>
        [MethodImpl(MethodImplOptions.Synchronized)]
        public void ProcessJobs(List<UnitOfWork> jobs, int maxNumThreads, bool redirectOutput, bool checkExitCodes, bool throwExceptionOnBadExitCode)
        {
            // sanity check: do we have any jobs?
            if (jobs.Count == 0) return;

            // we originally used TPL, but we encountered unexpected behavior from MaxDegreeOfParallelism
            // With MaxDegreeOfParallelism, it seemed like if it was set to 10 and each job used 50% CPU,
            // it would launch around 20 threads. ThreadPool was better, but also tried to be smart about
            // launching threads. Using standard threads and a semaphore yielded the desired behavior.

            // run our jobs
            _jobPool = new Semaphore(maxNumThreads, maxNumThreads);
            _doneEvent = new AutoResetEvent(false);
            _jobsRemaining = jobs.Count;
            bool cancelJobs = false;
            long maxVMem = 0;
            TimeSpan maxWallTime = new TimeSpan(0);
            for (int jobIndex = 0; jobIndex < jobs.Count; ++jobIndex)
            {
                _jobPool.WaitOne();

                // quit early if an exception was encountered in the thread
                if (cancelJobs)
                {
                    _doneEvent.Set();
                    break;
                }

                UnitOfWork job = jobs[jobIndex];
                Thread jobThread = new Thread(o => RunJob(job, redirectOutput, throwExceptionOnBadExitCode,
                                                          ref cancelJobs, ref maxVMem, ref maxWallTime));
                jobThread.Name = string.Format("Job thread {0}", Path.GetFileName(job.ExecutablePath));
                jobThread.Start();
            }

            _doneEvent.WaitOne();
            OnLog(String.Format("Max job memory (GB): {0:F1}", maxVMem / ((double)1024 * 1024 * 1024)));
            OnLog(String.Format("Longest job runtime (Hours): {0:F2}", maxWallTime.TotalHours));
            if (checkExitCodes) CheckExitCodes(jobs, throwExceptionOnBadExitCode);
        }

        /// <summary>
        /// processes a single job
        /// </summary>
        public void ProcessJob(UnitOfWork job, bool redirectOutput, bool checkExitCodes, bool throwExceptionOnBadExitCode)
        {
            ProcessJobs(new List<UnitOfWork> { job }, 1, redirectOutput, checkExitCodes, throwExceptionOnBadExitCode);
        }

        /// <summary>
        /// processes a single job
        /// </summary>
        private void RunJob(UnitOfWork job, bool redirectOutput, bool throwExceptionOnBadExitCode,
                            ref bool cancelJobs, ref long maxVMem, ref TimeSpan maxWallTime)
        {
            // replace the command-line tokens
            job.CommandLine = ConvertCommandLineTokens(job.CommandLine);
            string executableFilename = Path.GetFileName(job.ExecutablePath);

            OnLog(string.Format("Launch process: {0} {1}", job.ExecutablePath, job.CommandLine));

            // set the working directory to a drive we can write to (NOT a network drive)
            string workingDirectory = string.IsNullOrEmpty(job.WorkingDirectory) ? _alternativeWorkingDirectory : job.WorkingDirectory;

            foreach (string key in job.Environment.Keys)
            {
                job.CommandLine = job.CommandLine.Replace("%" + key + "%", job.Environment[key]);
            }

            ProcessStartInfo startInfo = new ProcessStartInfo(job.ExecutablePath, job.CommandLine)
            {
                WorkingDirectory = workingDirectory,
                CreateNoWindow = true,
                UseShellExecute = false
            };

            foreach (string key in job.EnvironmentVariables.Keys)
            {
                if (job.EnvironmentVariables[key] == null && startInfo.EnvironmentVariables.ContainsKey(key))
                    startInfo.EnvironmentVariables.Remove(key);
                else
                    startInfo.EnvironmentVariables[key] = job.EnvironmentVariables[key];
            }

            Process jobProcess = new Process { StartInfo = startInfo };

            // handle output redirection
            OutputRedirector outputDirector = null;
            OutputRedirector errorDirector = null;

            DeleteLog(job.ErrorLogPath);
            DeleteLog(job.OutputLogPath);
            //DeleteLog(job.CommandLogPath);

            if (redirectOutput)
            {
                startInfo.RedirectStandardError = true;
                startInfo.RedirectStandardOutput = true;

                outputDirector = new OutputRedirector(job.OutputLogPath);
                errorDirector = new OutputRedirector(job.ErrorLogPath);

                jobProcess.OutputDataReceived += outputDirector.LineHandler;
                jobProcess.ErrorDataReceived += errorDirector.LineHandler;
            }

            lock (_commandLogLock)
            {
                if (!string.IsNullOrEmpty(job.CommandLogPath))
                    using (var sw = new StreamWriter(job.CommandLogPath, true))
                    {
                        sw.WriteLine(job.LoggingStub + "\t" + job.ExecutablePath + " " + job.CommandLine);
                    }
            }

            // start the process
            bool stopMonitorThread = false;
            Thread monitorThread = null;
            try
            {
                jobProcess.Start();

                try
                {
                    if (Utilities.RunProcessesWithLowPriority) jobProcess.PriorityClass = ProcessPriorityClass.Idle;
                }
                catch (Exception)
                {
                    // the process may have completed so quickly we were unable to set the priority on it
                }

                // handle output redirection
                if (redirectOutput)
                {
                    jobProcess.BeginOutputReadLine();
                    jobProcess.BeginErrorReadLine();
                }

                // start the process monitor thread
                // todo monitoring could probably be done on this thread instead, slightly simplifying the code
                long procVMem = 0; //Max. mem used by this process (bytes)
                monitorThread = new Thread(o => MonitorProcess(jobProcess, ref stopMonitorThread, ref procVMem));
                monitorThread.Name = string.Format("Monitor thread {0}", executableFilename);
                monitorThread.Start();

                // wait for things to finish
                jobProcess.WaitForExit();
                lock (_rUsageLock)
                {
                    TimeSpan? wTime = GetProcessWallTime(jobProcess);
                    if (wTime.HasValue && wTime.Value > maxWallTime)
                        maxWallTime = wTime.Value;
                    if (procVMem > maxVMem) maxVMem = procVMem;
                }
                job.ExitCode = GetProcessExitCode(jobProcess);
            }
            catch (Exception e)
            {
                if (!throwExceptionOnBadExitCode) OnError(string.Format("Unable to start the process ({0}): {1}", executableFilename, e.Message));
                job.ThrownException = e;
                job.ExitCode = -1;
            }

            // cleanup
            stopMonitorThread = true;
            jobProcess.Close();
            jobProcess.Dispose();

            if (redirectOutput)
            {
                outputDirector.Dispose();
                errorDirector.Dispose();
                // Cleanup: Tidy up any dummy files!
                CleanDummyLogFile(job.OutputLogPath);
                CleanDummyLogFile(job.ErrorLogPath);
            }

            // prevent new jobs from running if we're going to throw an exception
            if (throwExceptionOnBadExitCode && (job.ExitCode != 0)) cancelJobs = true;

            if (job.ExitCode != 0) OnLog(string.Format("Process exited with code {2}: {0} {1}", job.ExecutablePath, job.CommandLine, job.ExitCode));

            // update the thread status
            if (Interlocked.Decrement(ref _jobsRemaining) <= 0) _doneEvent.Set();
            _jobPool.Release();
        }

        private void CleanDummyLogFile(string outputLogPath)
        {
            try
            {
                if (File.Exists(outputLogPath))
                {
                    FileInfo info = new FileInfo(outputLogPath);
                    if (info.Length == 0) // Ignore a file that consists only of \n or \r\n
                    {
                        File.Delete(outputLogPath);
                    }
                }
            }
            catch
            {
                ;
            }
        }

    }

    public class UnitOfWork
    {
        #region standard
        // input
        public string CommandLine;
        public string ExecutablePath;
        public string LoggingFolder;
        public string LoggingStub;

        // output
        public int ExitCode = -1; //Default to failure
        public string ErrorLogPath
        {
            get
            {
                if (string.IsNullOrEmpty(LoggingFolder) || string.IsNullOrEmpty(LoggingStub)) return null;
                return Path.Combine(LoggingFolder, LoggingStub) + ".stderr.txt";
            }
        }
        public string OutputLogPath
        {
            get
            {
                if (string.IsNullOrEmpty(LoggingFolder) || string.IsNullOrEmpty(LoggingStub)) return null;
                return Path.Combine(LoggingFolder, LoggingStub) + ".stdout.txt";
            }
        }

        public string CommandLogPath
        {
            get
            {
                if (string.IsNullOrEmpty(LoggingFolder)) return null;
                return Path.Combine(LoggingFolder, "commands.tsv");
            }
        }

        public Exception ThrownException = null;
        #endregion

        #region non-standard
        public string WorkingDirectory;
        public Dictionary<string, string> Environment = new Dictionary<string, string>(); //use this to customize a generic command line by replacing each appearance of %key% with a corresponding value
        public Dictionary<string, string> EnvironmentVariables = new Dictionary<string, string>(); //use this to set custom environment variables required by the process. to remove an environment variable set the value to null
        #endregion

        #region cluster-specific

        public string JobID;
        public double MaxMemoryUsageGB;

        #endregion
    }

    [Serializable]
    public class JobFailedException : Exception
    {
        public JobFailedException(string message)
            : base(message)
        {
        }
        protected JobFailedException(SerializationInfo info, StreamingContext ctxt) : base(info, ctxt)
        {
        }
    }

}
