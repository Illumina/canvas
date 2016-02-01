#region

using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Runtime.InteropServices;
using System.Security;
using System.Text;
using Illumina.SecondaryAnalysis;

#endregion

namespace Drmaa
{
    /// <summary>
    ///     This class handles the interoperability between C# and the drmaa C library
    ///     N.B. This class is not thread safe
    /// </summary>
    public class Interop : IDisposable
    {
        #region members

        private const string WrapperFilename = "SgeWrapper.exe";

        private const char SingleQuote = '\'';
        private const char DoubleQuote = '\"';
        private const char RecordSeparator = (char) 30;
        private const char EqualSign = '=';

        private readonly ErrorHandler onLog;
        private readonly ErrorHandler onError;
        private readonly string _parallelEnvironmentName;
        private readonly string _queueName;
        private readonly string _wrapperPath;
        private bool _isDisposed;
        private bool _isInitialized;

        #endregion

        #region DRMAA constants
        #pragma warning disable 169

        // from drmaa.h
        private const int DrmaaErrorStringBuffer = 1024;
        private const int DrmaaJobnameBuffer     = 1024;
        private const int DrmaaAttrBuffer        = 1024;

        private const int DrmaaTimeoutWaitForever = -1;
        private const int DrmaaTimeoutNoWait      = 0;

        private const int DrmaaErrnoSuccess        = 0;
        private const int DrmaaErrnoNoMoreElements = 25;

        private const string DrmaaVArgv  = "drmaa_v_argv";
        private const string DrmaaVEnv   = "drmaa_v_env";
        private const string DrmaaVEmail = "drmaa_v_email";

        private const string DrmaaRemoteCommand       = "drmaa_remote_command";
        private const string DrmaaJsState             = "drmaa_js_state";
        private const string DrmaaWd                  = "drmaa_wd";
        private const string DrmaaJobCategory         = "drmaa_job_category";
        private const string DrmaaNativeSpecification = "drmaa_native_specification";
        private const string DrmaaBlockEmail          = "drmaa_block_email";
        private const string DrmaaStartTime           = "drmaa_start_time";
        private const string DrmaaJobName             = "drmaa_job_name";
        private const string DrmaaInputPath           = "drmaa_input_path";
        private const string DrmaaOutputPath          = "drmaa_output_path";
        private const string DrmaaErrorPath           = "drmaa_error_path";
        private const string DrmaaJoinFiles           = "drmaa_join_files";
        private const string DrmaaTransferFiles       = "drmaa_transfer_files";
        private const string DrmaaDeadlineTime        = "drmaa_deadline_time";
        private const string DrmaaWctHlimit           = "drmaa_wct_hlimit";
        private const string DrmaaWctSlimit           = "drmaa_wct_slimit";
        private const string DrmaaDurationHlimit      = "drmaa_duration_hlimit";
        private const string DrmaaDurationSlimit      = "drmaa_duration_slimit";

        #pragma warning restore 169
        #endregion

        /// <summary>
        /// Wrapper for SGE job submissions
        /// </summary>
        /// <param name="queueName"> Queue to request (-q ...). Can be null to not request a specific queue</param>
        /// <param name="parallelEnvironmentName">Used to request slots (-pe ...)</param>
        /// <param name="log"></param>
        /// <param name="error"></param>
        public Interop(string queueName, string parallelEnvironmentName, ErrorHandler log, ErrorHandler error)
        {
            _queueName = queueName;
            _parallelEnvironmentName = parallelEnvironmentName;
            onLog = (log == null) ? Console.WriteLine : log;
            onError = (error == null) ? Console.WriteLine : error;

            // define the path to SgeWrapper.exe
            string executableFolder = Utilities.GetAssemblyFolder(typeof (Interop));
            _wrapperPath            = Path.Combine(executableFolder, WrapperFilename);

            if (!File.Exists(_wrapperPath))
            {
                throw new ApplicationException(string.Format("ERROR: Could not find the {0} executable ({1}).", WrapperFilename, _wrapperPath));
            }

            Initialize();
        }

        /// <summary>
        ///     starts a list of job on the cluster using SGE and then waits for the jobs to finish
        /// </summary>
        public void RunJobs(List<UnitOfWork> jobs, TaskResourceRequirements taskRequirements, bool throwExceptionOnError)
        {
            bool failed = false; // Did any job fail
            // sanity check: make sure we're initialized
            if (!_isInitialized)
            {
                throw new ApplicationException("Tried to run a job when the DRMAA library wasn't initialized yet.");
            }

            // sanity check: make sure we have jobs to execute
            if (jobs.Count == 0)
            {
                onLog("Skipping job submission to the cluster because no jobs were specified.");
                return;
            }

            // initialize
            IntPtr jobTemplate = AllocateJobTemplate();

            // submit all of the jobs
            foreach (UnitOfWork unitOfWork in jobs)
            {
				if (!string.IsNullOrEmpty(unitOfWork.CommandLogPath))
					using (var sw = new StreamWriter(unitOfWork.CommandLogPath, true))
					{
						sw.WriteLine(unitOfWork.LoggingStub + "\t" + unitOfWork.ExecutablePath + " " + unitOfWork.CommandLine);
					}
                // create the native specification
                string nativeSpecification = 
                    $"-V -pe {_parallelEnvironmentName} {taskRequirements.GetMinSlots()}-{taskRequirements.GetMaxSlots()}";
                if (!String.IsNullOrEmpty(_queueName)) nativeSpecification += $" -q {_queueName}";
                // set the job attributes
                SetAttribute(jobTemplate, DrmaaRemoteCommand, Utilities.GetMonoPath());
                SetAttribute(jobTemplate, DrmaaJobName, unitOfWork.LoggingStub);
                SetAttribute(jobTemplate, DrmaaOutputPath, ":" + unitOfWork.OutputLogPath);
                SetAttribute(jobTemplate, DrmaaErrorPath, ":" + unitOfWork.ErrorLogPath);
                SetAttribute(jobTemplate, DrmaaNativeSpecification, nativeSpecification);

                // handle the command-line arguments
                List<string> envList = new List<string>
                {
                    string.Format("ISIS_DRMAA_PATH={0}", unitOfWork.ExecutablePath),
                    string.Format("ISIS_DRMAA_ARGUMENTS={0}", unitOfWork.CommandLine.Replace(EqualSign, RecordSeparator)),
                    string.Format("ISIS_CLUSTER_CORES_PER_SLOT={0}", TaskResourceRequirements.CoresPerSGESlot),
                    string.Format("ISIS_CLUSTER_GB_PER_SLOT={0}", TaskResourceRequirements.GigabytesPerSGESlot),
                    string.Format("ISIS_MIN_CORES={0}", taskRequirements.MinCores),
                    string.Format("ISIS_MAX_CORES={0}", taskRequirements.MaxCores),
                    string.Format("ISIS_MIN_GB={0}", taskRequirements.MinMemoryGB),
                    string.Format("ISIS_EXTRA_GB_PER_CORE={0}", taskRequirements.MemoryGBPerExtraCore),
                    string.Format("ISIS_WORKING_DIRECTORY={0}", unitOfWork.WorkingDirectory),
                    string.Format("ISIS_SGE_COMMAND={0}", nativeSpecification),
                };

                string[] environmentVariables = envList.ToArray();
                SetAttribute(jobTemplate, DrmaaVEnv, environmentVariables);

                string[] arguments = new[] {_wrapperPath};
                SetAttribute(jobTemplate, DrmaaVArgv, arguments);

                // run the job
                unitOfWork.JobID = ExecuteJob(jobTemplate);
                if (unitOfWork.JobID == null) failed = true;
                onLog(string.Format("Submit: {0} {1}", unitOfWork.ExecutablePath, unitOfWork.CommandLine));
            }

            string logStr = $"{jobs.Count} jobs have been submitted";
            if (!String.IsNullOrEmpty(_queueName)) logStr += $" to queue {_queueName}";
            logStr += String.Format(". SGE slots requested: {0}:{1}-{2}", _parallelEnvironmentName, 
                    taskRequirements.GetMinSlots(), taskRequirements.GetMaxSlots());
            onLog(logStr);
            // synchronize all jobs
            if (!WaitForJobs(jobs)) failed = true;

            onLog("All cluster jobs have now completed.");

            // grab the resource usage for each job
            double maxMem = 0; // Memory usage of largest job
            TimeSpan maxTime = new TimeSpan(0); // Runtime of longest job
            foreach (UnitOfWork unitOfWork in jobs)
            {
                IntPtr resourceUsage;
                if (!WaitForJob(unitOfWork, out resourceUsage)) failed = true; //todo Cancel other SGE jobs here if requested
                // save the resource usage
                DrmaaResourceUsage rUsage = GetResourceUsage(resourceUsage);
                SaveResourceUsage(rUsage, unitOfWork);
                maxMem = Math.Max(maxMem, rUsage.MaximumVirtualMemory);
                if (rUsage.WallTime > maxTime) maxTime = rUsage.WallTime;
                // cleanup
                ReleaseAttributeValues(resourceUsage);
            }
            onLog(String.Format("Max job memory (GB): {0:F1}", maxMem / ((double)1024 * 1024 * 1024)));
            onLog(String.Format("Longest job runtime (Hours): {0:F2}", maxTime.TotalHours));
            if (failed)
            {
                onError("At least one SGE job failed");
                if (throwExceptionOnError) throw new JobFailedException("SGE job failed");
            }
            // cleanup
            DeleteJobTemplate(jobTemplate);
        }

        #region IDisposable

        // Implement IDisposable. 
        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        // Dispose(bool disposing) executes in two distinct scenarios. 
        // If disposing equals true, the method has been called directly 
        // or indirectly by a user's code. Managed and unmanaged resources 
        // can be disposed. 
        // If disposing equals false, the method has been called by the 
        // runtime from inside the finalizer and you should not reference 
        // other objects. Only unmanaged resources can be disposed. 
        protected virtual void Dispose(bool disposing)
        {
            // Check to see if Dispose has already been called. 
            if (!_isDisposed)
            {
                // Dispose managed resources.
                if (disposing)
                {
                    //component.Dispose();
                }

                Shutdown();

                //// Call the appropriate methods to clean up 
                //// unmanaged resources here. 
                //// If disposing is false, 
                //// only the following code is executed.
                //CloseHandle(handle);
                //handle = IntPtr.Zero;

                _isDisposed = true;
            }
        }

        #endregion

        #region Helper Functions

        /// <summary>
        ///     allocates a job template
        /// </summary>
        private IntPtr AllocateJobTemplate()
        {
            IntPtr jobTemplate = IntPtr.Zero;
            byte[] errorMessageBuffer = new byte[DrmaaErrorStringBuffer + 1];

            int errnum = SafeNativeMethods.drmaa_allocate_job_template(ref jobTemplate, errorMessageBuffer,
                                                                       DrmaaErrorStringBuffer);

            if (errnum != DrmaaErrnoSuccess)
            {
                string error = ConvertByteArrayToString(errorMessageBuffer);
                throw new ApplicationException(string.Format("Could not create job template: {0}", error));
            }

            return jobTemplate;
        }

        /// <summary>
        ///     returns a string given a null-terminated byte array
        /// </summary>
        private static string ConvertByteArrayToString(byte[] byteArray)
        {
            // retrieve the location of the null
            int nullIndex = Array.IndexOf<byte>(byteArray, 0);

            if (nullIndex == -1)
            {
                throw new ApplicationException("Unable to find the null-termination in the byte array.");
            }

            // convert the byte array into a char array
            char[] s = new char[nullIndex];
            for (int i = 0; i < nullIndex; ++i) s[i] = (char) byteArray[i];

            // convert the char array to a string
            return new string(s);
        }

        /// <summary>
        ///     frees up the resources used by this job template
        /// </summary>
        private void DeleteJobTemplate(IntPtr jobTemplate)
        {
            byte[] errorMessageBuffer = new byte[DrmaaErrorStringBuffer + 1];

            int errnum = SafeNativeMethods.drmaa_delete_job_template(jobTemplate, errorMessageBuffer,
                                                                     DrmaaErrorStringBuffer);

            if (errnum != DrmaaErrnoSuccess)
            {
                string error = ConvertByteArrayToString(errorMessageBuffer);
                throw new ApplicationException(string.Format("Could not delete job template: {0}", error));
            }
        }

        /// <summary>
        ///     executes the job defined in the job template
        /// </summary>
        /// <returns>the job ID</returns>
        private string ExecuteJob(IntPtr jobTemplate)
        {
            byte[] jobNameBuffer      = new byte[DrmaaJobnameBuffer + 1];
            byte[] errorMessageBuffer = new byte[DrmaaErrorStringBuffer + 1];

            int errNum = SafeNativeMethods.drmaa_run_job(jobNameBuffer, DrmaaJobnameBuffer, jobTemplate,
                                                         errorMessageBuffer, DrmaaErrorStringBuffer);

            if (errNum != DrmaaErrnoSuccess)
            {
                string error = ConvertByteArrayToString(errorMessageBuffer);
                onError(string.Format("Could not run the job: {0}", error));
                return null;
            }

            return ConvertByteArrayToString(jobNameBuffer);
        }

        /// <summary>
        ///     return a DateTime given a double in string form
        /// </summary>
        private static DateTime GetRusageDateTime(Dictionary<string, string> resourceData, string key)
        {
            DateTime dt = new DateTime(1970, 1, 1, 0, 0, 0, 0);
            dt = dt.AddSeconds(GetRusageDouble(resourceData, key)).ToLocalTime();
            return dt;
        }

        /// <summary>
        ///     return an integer given a double in string form
        /// </summary>
        private static double GetRusageDouble(Dictionary<string, string> resourceData, string key)
        {
            string value;
            if (!resourceData.TryGetValue(key, out value)) return default(double);

            double d;
            if (!double.TryParse(value, out d))
            {
                throw new ApplicationException(string.Format("Unable to convert the string ({0}) to a double.", value));
            }

            return d;
        }

        /// <summary>
        ///     return an integer given a double in string form
        /// </summary>
        private static int GetRusageInt(Dictionary<string, string> resourceData, string key)
        {
            return (int) GetRusageDouble(resourceData, key);
        }

        /// <summary>
        ///     return a TimeSpan given a double in string form
        /// </summary>
        private static TimeSpan GetRusageTimeSpan(Dictionary<string, string> resourceData, string key)
        {
            return TimeSpan.FromSeconds(GetRusageDouble(resourceData, key));
        }

        /// <summary>
        ///     populates a resource usage data structure with the resource usage statistics from this job
        /// </summary>
        private DrmaaResourceUsage GetResourceUsage(IntPtr resourceUsage)
        {
            // initialize
            byte[] attributeBuffer                  = new byte[DrmaaAttrBuffer + 1];
            Dictionary<string, string> resourceData = new Dictionary<string, string>();

            // put all of the attributes in a dictionary
            while (SafeNativeMethods.drmaa_get_next_attr_value(resourceUsage, attributeBuffer, DrmaaAttrBuffer) !=
                   DrmaaErrnoNoMoreElements)
            {
                string[] cols = ConvertByteArrayToString(attributeBuffer).Split('=');
                resourceData[cols[0]] = (cols.Length == 2 ? cols[1] : string.Empty);
            }

            // populate the resource usage data structure
            DrmaaResourceUsage rusage = new DrmaaResourceUsage
                {
                    CpuTime              = GetRusageTimeSpan(resourceData, "cpu"),
                    UserTime             = GetRusageTimeSpan(resourceData, "ru_utime"),
                    IoWaitTime           = GetRusageTimeSpan(resourceData, "iow"),
                    WallTime             = GetRusageTimeSpan(resourceData, "ru_wallclock"),
                    SubmissionTime       = GetRusageDateTime(resourceData, "submission_time"),
                    StartTime            = GetRusageDateTime(resourceData, "start_time"),
                    EndTime              = GetRusageDateTime(resourceData, "end_time"),
                    IoDataTransferredGB  = GetRusageDouble(resourceData, "io"),
                    NumPageFaults        = GetRusageInt(resourceData, "ru_majflt"),
                    NumPageReclaims      = GetRusageInt(resourceData, "ru_minflt"),
                    VirtualMemory        = GetRusageDouble(resourceData, "vmem"),
                    MaximumVirtualMemory = GetRusageDouble(resourceData, "maxvmem"),
                    Priority             = GetRusageDouble(resourceData, "priority"),
                    ExitStatus           = GetRusageInt(resourceData, "exit_status"),
                    Signal               = GetRusageInt(resourceData, "signal")
                };

            return rusage;
        }

        /// <summary>
        ///     initializes the DRMAA library functions
        /// </summary>
        private void Initialize()
        {
            byte[] errorMessageBuffer = new byte[DrmaaErrorStringBuffer + 1];

            int errNum = SafeNativeMethods.drmaa_init(IntPtr.Zero, errorMessageBuffer, DrmaaErrorStringBuffer);

            if (errNum != DrmaaErrnoSuccess)
            {
                string error = ConvertByteArrayToString(errorMessageBuffer);
                throw new ApplicationException(string.Format("Could not initialize the DRMAA library: {0}", error));
            }

            _isInitialized = true;
        }

        /// <summary>
        ///     frees up the resources used by the resource usage attributes
        /// </summary>
        private void ReleaseAttributeValues(IntPtr resourceUsage)
        {
            SafeNativeMethods.drmaa_release_attr_values(resourceUsage);
        }

        /// <summary>
        ///     saves the output from the resource usage object to an output file
        /// </summary>
        private void SaveResourceUsage(DrmaaResourceUsage rUsage, UnitOfWork unitOfWork)
        {
            // derive the resource usage filename
            string resourcesPath = string.Format("{0}.cluster", Path.Combine(unitOfWork.LoggingFolder, unitOfWork.LoggingStub));

            // save the resource usage to the output file
            try
            {
                using (StreamWriter writer = new StreamWriter(resourcesPath))
                {
                    writer.WriteLine(rUsage);
                }
            }
            catch (Exception)
            {
                // this is a non-fatal error
                Console.WriteLine("ERROR: Unable to write the resource usage file ({0}).", resourcesPath);
            }
        }

        /// <summary>
        ///     sets the key/value data in the job template
        /// </summary>
        private void SetAttribute(IntPtr jobTemplate, string key, string value)
        {
            byte[] errorMessageBuffer = new byte[DrmaaErrorStringBuffer + 1];

            ASCIIEncoding encoding = new ASCIIEncoding();
            byte[] keyBuffer       = encoding.GetBytes(key + char.MinValue);
            byte[] valueBuffer     = encoding.GetBytes(value + char.MinValue);

            int errNum = SafeNativeMethods.drmaa_set_attribute(jobTemplate, keyBuffer, valueBuffer, errorMessageBuffer, DrmaaErrorStringBuffer);

            if (errNum != DrmaaErrnoSuccess)
            {
                string error = ConvertByteArrayToString(errorMessageBuffer);
                throw new ApplicationException(string.Format("Could not set attribute ({0}) = ({1}): {2}", key, value, error));
            }
        }

        /// <summary>
        ///     sets the key/value data in the job template
        /// </summary>
        private void SetAttribute(IntPtr jobTemplate, string key, string[] value)
        {
            byte[] errorMessageBuffer = new byte[DrmaaErrorStringBuffer + 1];
            ASCIIEncoding encoding    = new ASCIIEncoding();
            byte[] keyBuffer          = encoding.GetBytes(key + char.MinValue);

            // the expectation here is that the array terminates with a null character
            string[] valueWithNull = new string[value.Length + 1];
            value.CopyTo(valueWithNull, 0);
            valueWithNull[value.Length] = null;

            int errNum = SafeNativeMethods.drmaa_set_vector_attribute(jobTemplate, keyBuffer, valueWithNull, errorMessageBuffer, DrmaaErrorStringBuffer);

            if (errNum != DrmaaErrnoSuccess)
            {
                string error = ConvertByteArrayToString(errorMessageBuffer);
                throw new ApplicationException(string.Format("Could not set vector attribute ({0}): {1}", key, error));
            }
        }

        /// <summary>
        ///     shuts down the DRMAA library functions
        /// </summary>
        private void Shutdown()
        {
            byte[] errorMessageBuffer = new byte[DrmaaErrorStringBuffer + 1];

            int errnum = SafeNativeMethods.drmaa_exit(errorMessageBuffer, DrmaaErrorStringBuffer);

            if (errnum != DrmaaErrnoSuccess)
            {
                string error = ConvertByteArrayToString(errorMessageBuffer);
                throw new ApplicationException(string.Format("Could not shut down the DRMAA library: {0}", error));
            }

            _isInitialized = false;
        }

        /// <summary>
        ///     waits for the specified job to finish
        /// </summary>
        private bool WaitForJob(UnitOfWork unitOfWork, out IntPtr resourceUsage)
        {
            byte[] errorMessageBuffer = new byte[DrmaaErrorStringBuffer + 1];
            int status                = -1;
            resourceUsage             = IntPtr.Zero;

            // convert the job ID
            ASCIIEncoding encoding  = new ASCIIEncoding();
            byte[] jobNameBuffer    = encoding.GetBytes(unitOfWork.JobID + char.MinValue);
            byte[] jobNameOutBuffer = new byte[DrmaaJobnameBuffer + 1];

            int errnum = SafeNativeMethods.drmaa_wait(jobNameBuffer, jobNameOutBuffer, DrmaaJobnameBuffer, ref status,
                                                      DrmaaTimeoutWaitForever, ref resourceUsage, errorMessageBuffer,
                                                      DrmaaErrorStringBuffer);

            if (errnum != DrmaaErrnoSuccess)
            {
                string error = ConvertByteArrayToString(errorMessageBuffer);
                onError(string.Format("Could not wait for job: {0}", error));
                return false;
            }

            // Use the DRMAA function call to resolve the value of status
            int isExited = 0;
            SafeNativeMethods.drmaa_wifexited(ref isExited, status, errorMessageBuffer, DrmaaErrorStringBuffer);

            // isExited evaluates to zero if the job did not terminate normally
            if (isExited == 0)
            {
                string error = ConvertByteArrayToString(errorMessageBuffer);
                if (error == null || error == "") error = "<MISSING SGE ERROR>";
                onError(string.Format("Job exited abnormally: {0}.\n** Failed command:\n{1} {2}\n** Stderr output (if any) is at: {3}\n",
                    error, unitOfWork.ExecutablePath, unitOfWork.CommandLine, unitOfWork.ErrorLogPath));
                return false;
            }

            // Job finished normally, but maybe exited with a non-zero exit code, 
            // we'll check for that now
            unitOfWork.ExitCode = 0;
            SafeNativeMethods.drmaa_wexitstatus(ref unitOfWork.ExitCode, status, errorMessageBuffer, DrmaaErrorStringBuffer);
            if (unitOfWork.ExitCode != 0)
            {
                string error = ConvertByteArrayToString(errorMessageBuffer);
                onError(string.Format("Job returned non-zero exit code {0}: {1}.\n** Failed command:\n{2} {3}\n** Stderr output (if any) is at: {4}\n",
                    unitOfWork.ExitCode, error, unitOfWork.ExecutablePath, unitOfWork.CommandLine, unitOfWork.ErrorLogPath));
                return false;
            }
            onLog(string.Format("Successfully finished command:\n{0} {1}", unitOfWork.ExecutablePath, unitOfWork.CommandLine));
            return true;
        }

        /// <summary>
        ///     waits for a number of jobs to finish
        /// </summary>
        private bool WaitForJobs(List<UnitOfWork> jobs)
        {
            byte[] errorMessageBuffer = new byte[DrmaaErrorStringBuffer + 1];

            // grab the jobs IDs
            // N.B. the expectation here is that the array terminates with a null character
            string[] jobIdsWithNull = new string[jobs.Count + 1];

            for (int jobIndex = 0; jobIndex < jobs.Count; ++jobIndex)
            {
                jobIdsWithNull[jobIndex] = jobs[jobIndex].JobID;
            }

            jobIdsWithNull[jobs.Count] = null;

            // wait for the jobs to finish
            int errnum = SafeNativeMethods.drmaa_synchronize(jobIdsWithNull, DrmaaTimeoutWaitForever, 0,
                                                             errorMessageBuffer, DrmaaErrorStringBuffer);

            if (errnum != DrmaaErrnoSuccess)
            {
                string error = ConvertByteArrayToString(errorMessageBuffer);
                onError(string.Format("Could not synchronize jobs: {0}", error));
                return false;
            }
            return true;
        }

        #endregion
    }

    [SuppressUnmanagedCodeSecurity]
    internal static class SafeNativeMethods
    {
        // ReSharper disable InconsistentNaming
        [DllImport("drmaa.dll", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int drmaa_init(IntPtr contact, byte[] errorDiagnosis, int errorDiagLen);

        [DllImport("drmaa.dll", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int drmaa_exit(byte[] errorDiagnosis, int errorDiagLen);

        [DllImport("drmaa.dll", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int drmaa_allocate_job_template(ref IntPtr jt, byte[] errorDiagnosis, int errorDiagLen);

        [DllImport("drmaa.dll", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int drmaa_set_attribute(IntPtr jt, byte[] name, byte[] value, byte[] errorDiagnosis,
                                                       int errorDiagLen);

        [DllImport("drmaa.dll", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int drmaa_run_job(byte[] jobId, int jobIdLen, IntPtr jt, byte[] errorDiagnosis,
                                                 int errorDiagLen);

        [DllImport("drmaa.dll", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int drmaa_wait(byte[] jobId, byte[] jobIdOut, int jobIdOutLen, ref int stat,
                                              long timeout, ref IntPtr rusage, byte[] errorDiagnosis,
                                              int errorDiagLen);

        [DllImport("drmaa.dll", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int drmaa_delete_job_template(IntPtr jt, byte[] errorDiagnosis, int errorDiagLen);

        [DllImport("drmaa.dll", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int drmaa_get_next_attr_value(IntPtr values, byte[] value, int valueLen);

        [DllImport("drmaa.dll", CallingConvention = CallingConvention.Cdecl)]
        internal static extern void drmaa_release_attr_values(IntPtr values);

        [DllImport("drmaa.dll", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
        internal static extern int drmaa_set_vector_attribute(IntPtr jt, byte[] name, string[] value,
                                                              byte[] errorDiagnosis, int errorDiagLen);

        [DllImport("drmaa.dll", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
        internal static extern int drmaa_synchronize(string[] jobIds, long timeout, int dispose, byte[] errorDiagnosis,
                                                     int errorDiagLen);

        [DllImport("drmaa.dll", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
        internal static extern int drmaa_wifexited(ref int exited, int stat, byte [] error_diagnosis, int error_diag_len);

        [DllImport("drmaa.dll", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
        internal static extern int drmaa_wexitstatus(ref int exit_status, int stat, byte [] error_diagnosis, int error_diag_len);
        // ReSharper restore InconsistentNaming
        
        //public static int drmaa_wifaborted(int *aborted, int stat, char *error_diagnosis, size_t error_diag_len);
        //public static int drmaa_wexitstatus(int *exit_status, int stat, char *error_diagnosis, size_t error_diag_len);
        //public static int drmaa_wifsignaled(int *signaled, int stat, char *error_diagnosis, size_t error_diag_len);
        //public static int drmaa_wtermsig(char *signal, size_t signal_len, int stat, char *error_diagnosis, size_t error_diag_len);
    }

    internal class DrmaaResourceUsage
    {
        #region members

        // time
        private const double NumBytesIn1GB = 1073741824.0;
        private const double NumBytesIn1MB = 1048576.0;
        private const double NumBytesIn1KB = 1024.0;
        public TimeSpan CpuTime;
        public DateTime EndTime;
        public int ExitStatus;

        // I/O
        public double IoDataTransferredGB;
        public TimeSpan IoWaitTime;
        public double MaximumVirtualMemory;
        public int NumPageFaults;
        public int NumPageReclaims;

        // memory
        public double Priority;
        public int Signal;
        public DateTime StartTime;
        public DateTime SubmissionTime;
        public TimeSpan UserTime;
        public double VirtualMemory;
        public TimeSpan WallTime;

        #endregion

        /// <summary>
        ///     Creates a string representation of the resource usage object
        /// </summary>
        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            sb.AppendLine("Resource Usage:");
            sb.AppendLine("===============");
            sb.AppendLine();

            sb.AppendLine("Time        HH:MM:SS");
            sb.AppendLine("--------------------");
            if (CpuTime    != default(TimeSpan)) sb.AppendFormat("CPU:        {0}\n", GetTimeSpan(CpuTime));
            if (UserTime   != default(TimeSpan)) sb.AppendFormat("User:       {0}\n", GetTimeSpan(UserTime));
            if (WallTime   != default(TimeSpan)) sb.AppendFormat("Wall:       {0}\n", GetTimeSpan(WallTime));
            if (IoWaitTime != default(TimeSpan)) sb.AppendFormat("I/O Wait:   {0}\n", GetTimeSpan(IoWaitTime));
            sb.AppendLine();

            if (SubmissionTime != default(DateTime)) sb.AppendFormat("Submission: {0}\n", SubmissionTime);
            if (StartTime      != default(DateTime)) sb.AppendFormat("Start:      {0}\n", StartTime);
            if (EndTime        != default(DateTime)) sb.AppendFormat("End:        {0}\n", EndTime);
            sb.AppendLine();

            sb.AppendLine("Memory:");
            sb.AppendLine("-------");
            // ReSharper disable CompareOfFloatsByEqualityOperator
            if (VirtualMemory        != default(double)) sb.AppendFormat("Virtual:         {0}\n", GetMemoryUsage(VirtualMemory));
            if (MaximumVirtualMemory != default(double)) sb.AppendFormat("Maximum Virtual: {0}\n", GetMemoryUsage(MaximumVirtualMemory));
            // ReSharper restore CompareOfFloatsByEqualityOperator
            sb.AppendLine();

            sb.AppendLine("Job Info:");
            sb.AppendLine("---------");
            sb.AppendFormat("Priority:    {0}\n", Priority);
            sb.AppendFormat("Exit Status: {0}\n", ExitStatus);
            sb.AppendFormat("Signal:      {0}\n", Signal);

            return sb.ToString();
        }

        /// <summary>
        ///     returns a human-readable representation of how much memory was used
        /// </summary>
        private static string GetMemoryUsage(double numBytes)
        {
            // try GB
            double numGB = numBytes/NumBytesIn1GB;
            if (numGB > 1.0) return string.Format("{0:0.0} GB", numGB);

            // try MB
            double numMB = numBytes/NumBytesIn1MB;
            if (numMB > 1.0) return string.Format("{0:0.0} MB", numMB);

            // try KB
            double numKB = numBytes/NumBytesIn1KB;
            if (numKB > 1.0) return string.Format("{0:0.0} KB", numKB);

            // settle for B
            return string.Format("{0:0.} B", numBytes);
        }

        /// <summary>
        ///     returns a string depicting a TimeSpan object
        /// </summary>
        private static string GetTimeSpan(TimeSpan ts)
        {
            if (ts.Days > 0)
                return string.Format("{0:D2}d {1:D2}:{2:D2}:{3:D2}", ts.Days, ts.Hours, ts.Minutes, ts.Seconds);
            return string.Format("{0:D2}:{1:D2}:{2:D2}", ts.Hours, ts.Minutes, ts.Seconds);
        }
    }
}