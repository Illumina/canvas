using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading;
using System.Runtime.Serialization;
using System.Xml.Serialization;
using System.Runtime.Serialization.Formatters.Binary;
using Illumina.Common;
using Illumina.Zlib;
using ILMNcommon.Common;

namespace Isas.Shared
{
    public static class CommandLineTokens
    {
        public const string NumProcessors = "[NUM_PROCESSORS]";
    }

    public class OutputRedirector : IDisposable
    {
        #region Members

        private readonly string _filePath;
        private StreamWriter _writer;

        #endregion

        public OutputRedirector(string filePath, bool autoFlush = true)
        {
            _filePath = filePath;
            _writer = new StreamWriter(filePath);
            _writer.NewLine = "\n";
            _writer.AutoFlush = autoFlush;
        }

        public void Dispose()
        {
            if (_writer != null)
            {
                _writer.Close();
                _writer = null;
            }
        }

        public void LineHandler(object sendingProcess, DataReceivedEventArgs arguments)
        {
            if (_writer == null || arguments == null) return; // Just-in-case sanity check
            if (arguments.Data == null) return;
            try
            {
                _writer.WriteLine(arguments.Data);
            }
            catch (ObjectDisposedException)
            {
                // in case the TextWriter is closed already
                return;
            }
            catch (Exception ex)
            {
                Console.WriteLine(string.Format("Error writing to file: {0}\n  {1}", _filePath, ex));
            }
        }
    }

    public delegate void ErrorHandler(string message);

    public static class Utilities
    {
        public static bool RunProcessesWithLowPriority = true;
        private static readonly char[] CSVEscapeChars = new[] { ',', '"' };

        /// <summary>
        /// Return an iso-8601-formatted date.
        /// </summary>
        public static string FormatDate(DateTime date)
        {
            return date.ToString("yyyy-MM-dd", System.Globalization.CultureInfo.InvariantCulture);
        }

        /// <summary>
        /// Return an iso-8601-formatted date+time 
        /// </summary>
        public static string FormatDateTime(DateTime date)
        {
            return date.ToString("yyyy-MM-ddTHH:mm:ss", System.Globalization.CultureInfo.InvariantCulture);
        }

        /// <summary>
        /// Return an iso-8601-formatted date+time in UTC time (including trailing Z to indicate UTC time zone)
        /// </summary>
        public static string FormatUTCNow()
        {
            return DateTime.UtcNow.ToString("yyyy-MM-ddTHH:mm:ssZ", System.Globalization.CultureInfo.InvariantCulture);
        }

        //a hashset to check if a variable is of numeric type.
        // example:
        //  bool isANumber = NumericTypes.Contains(classInstance.GetType());
        public static HashSet<Type> NumericTypes = new HashSet<Type>
        {
            typeof(decimal), typeof(byte), typeof(sbyte),
            typeof(short), typeof(ushort),typeof(double),
            typeof(float),typeof(int), typeof(uint),typeof(Int64), typeof(UInt64)
        };

        public static string EscapeCommas(string theString)
        {
            if (theString == null) return "";
            if (theString.Contains("\"")) theString = theString.Replace("\"", "\"\"");
            if (theString.IndexOfAny(CSVEscapeChars) > -1) theString = "\"" + theString + "\"";
            return theString;
        }

        /// <summary>
        /// When calling a Windows command line and putting arguments between double quotation marks,
        ///   one can run in the issues of escaping the end quotation mark instead of encapsulating it.
        /// For example: 
        ///   pProcess.StartInfo.Arguments = String.Format("-a \"{0}\" -b \"{0}\"", arg1, arg2);
        ///   would cause issues if either arg1 or arg2 ended in a backslash (e.g. if it were a folder path)
        /// See CommandLineToArgvW for more details. 
        /// </summary>
        public static string EscapeTrailingBackslashes(string theString)
        {
            // According to MSDN website: 
            //     Backslashes are interpreted literally, unless they immediately precede a double quotation mark.
            //     See CommandLineToArgvW for more details.
            int numEndBackslashes = theString.Length - theString.TrimEnd('\\').Length;
            if (numEndBackslashes > 0)
            {
                theString += new String('\\', numEndBackslashes);
            }
            return theString;
        }

        public static bool IsThisMono()
        {
            return Type.GetType("Mono.Runtime") != null;
        }

        public enum TabixFileType
        {
            Bed,
            Vcf
        }

        public static string Preset(this TabixFileType type)
        {
            switch (type)
            {
                case TabixFileType.Bed:
                    return "bed";
                case TabixFileType.Vcf:
                    return "vcf";
                default:
                    throw new ArgumentException($"Unsupported tabix file type {type}");
            }
        }

        public static void BuildBedTabixIndex(string bedPath, string logFolder)
        {
            BuildBedTabixIndex(bedPath.ToSingleItemEnumerable(), logFolder);
        }

        public static void BuildBedTabixIndex(IEnumerable<string> bedPaths, string logFolder)
        {
            BuildTabixIndex(bedPaths, logFolder, TabixFileType.Bed);
        }

        public static void BuildVcfTabixIndex(string vcfFile, string logFolder)
        {
            BuildVcfTabixIndex(vcfFile.ToSingleItemEnumerable(), logFolder);
        }

        public static void BuildVcfTabixIndex(IEnumerable<string> vcfFiles, string logFolder)
        {
            BuildTabixIndex(vcfFiles, logFolder, TabixFileType.Vcf);
        }

        private static void BuildTabixIndex(IEnumerable<string> files, string logFolder, TabixFileType type)
        {
            // Look for tabix in worker folder, then in $PATH:
            string tabixExecutable;
            try
            {
                tabixExecutable = GetExecutablePath("tabix", null);
            }
            catch
            {
                string exeName = IsThisMono() ? "tabix" : "tabix.exe";
                tabixExecutable = CheckForExecutableInPath(exeName);
            }

            List<UnitOfWork> jobs = new List<UnitOfWork>();
            foreach (string fn in files)
            {
                if (!File.Exists(fn)) { continue; }
                string cmd = $" -p {type.Preset()} -f {fn.WrapWithShellQuote()}";
                jobs.Add(new UnitOfWork
                {
                    ExecutablePath = tabixExecutable,
                    LoggingFolder = logFolder,
                    LoggingStub = "tabix-" + Path.GetFileName(fn), // Note: Assume we only need the most recent logging if we end up tabix-ing the same file 2+ times
                    CommandLine = cmd
                });
            }
            Benchmark variantIndexingBenchmark = new Benchmark();
            DoWorkParallel(IsasConfiguration.GetConfiguration(), Console.WriteLine, Console.Error.WriteLine, jobs, new TaskResourceRequirements(1, 4));
            Console.WriteLine($"Total time to create tabix indexes for {files.Count()} files: {variantIndexingBenchmark.GetElapsedTime()}");
        }

        public static string GetJavaExecutable()
        {
            string javaExecutableName = "java";
            if (!IsThisMono()) javaExecutableName += ".exe";

            string executableFolder = GetAssemblyFolder(typeof(Utilities));
            string rootFolder = Path.GetDirectoryName(Path.GetDirectoryName(executableFolder));
            string javaPath = Path.Combine(rootFolder, "Java");
            string javaBinPath = Path.Combine(javaPath, "bin");
            string javaExecutablePath = Path.Combine(javaBinPath, javaExecutableName);

            if (!File.Exists(javaExecutablePath))
            {
                Console.WriteLine("Could not find the java executable ({0}). Searching in PATH for {1}.", javaExecutablePath, javaExecutableName);
                return FindExePath(javaExecutableName);
            }

            return javaExecutablePath;
        }

        /// <summary>
        ///     Convenience method for serializing an object to XML
        /// </summary>
        public static void XMLSerialize(object theThing, string filePath)
        {
            XmlSerializer xs = new XmlSerializer(theThing.GetType());
            using (FileStream fs = new FileStream(filePath, FileMode.Create))
            {
                xs.Serialize(fs, theThing);
                fs.Flush(); //flush it
                fs.Close(); //close it (to avoid Sharing violation in moving this file right afterwards)
            }
        }

        /// <summary>
        /// Quick-and-dirty bash shell call
        /// </summary>
        public static void ExecuteShellCommand(string cmd)
        {
            ExecuteCommand("bash", String.Format("-c '{0}'", cmd));
        }

        public static UnitOfWork CreateOneUnitOfWork(string executablePath = null, string commandLineArguments = "",
            StringBuilder commandBuilder = null,
            string loggingStub = "", string loggingFolder = "")
        {
            if (string.IsNullOrEmpty(executablePath))
            {
                throw new Exception(string.Format("executable {0} not found.", executablePath));
            }
            if (commandBuilder != null)
            {
                commandLineArguments = commandLineArguments + commandBuilder.ToString();
            }
            var job = new UnitOfWork
            {
                ExecutablePath = executablePath,
                LoggingFolder = loggingFolder,
                LoggingStub = loggingStub,
                CommandLine = commandLineArguments
            };

            return job;
        }

        /// <summary>
        ///  Try several times to delete a given directory tree (sometimes needed on NFS file-systems).
        /// </summary>
        public static bool SafeDeleteRetry(string path, int retries = 3, int sleep = 10000)
        {
            if (!Directory.Exists(path)) return false;
            bool deleted = false;
            for (int retryCount = 0; retryCount < retries; retryCount += 1)
            {
                try
                {
                    Directory.Delete(path, true);
                    deleted = true;
                    break;
                }
                catch (IOException)
                {
                }
                Thread.Sleep(sleep);
            }
            return deleted;
        }

        /// <summary>
        /// Quick-and-dirty single command-call without stdout/stderr
        /// </summary>
        public static void ExecuteCommand(string executable, string arguments, ErrorHandler onLog = null, ErrorHandler onError = null)
        {
            if (onLog == null)
                onLog = s => Console.WriteLine(s);
            if (onError == null)
                onError = s => Console.Error.WriteLine(s);
            UnitOfWork job = new UnitOfWork
            {
                ExecutablePath = executable,
                CommandLine = arguments,
            };
            DoWorkParallelThreads(IsasConfiguration.GetConfiguration(), onLog, onError, new List<UnitOfWork> { job }, 1, true, false, false);
        }

        // BinaryFormatter might not be thread safe so we lock on this object
        private static Object _binarySerializationLock = new Object();

        /// <summary>
        ///     Convenience method for serializing an object
        ///     BinaryFormatter allows you to serialize private member data unlike XmlSerializer
        /// </summary>
        public static void BinarySerialize(object TheThing, string FilePath)
        {
            lock (_binarySerializationLock)
            {
                BinaryFormatter xs = new BinaryFormatter();
                // allow shared read/write access for deserialization in case file lock is not released even though this serialization process has exited
                using (FileStream fs = new FileStream(FilePath, FileMode.Create, FileAccess.Write, FileShare.ReadWrite))
                {
                    xs.Serialize(fs, TheThing);
                }
            }
        }

        // hack for missing assembly exception during binary deserialization in .NET
        // maybe .NET does not like assemblies being copied around? Mono doesn't have a problem finding the assembly
        public class MySerializationBinder : SerializationBinder
        {
            private Type T;
            public MySerializationBinder(Type t)
            {
                T = t;
            }
            public override Type BindToType(string assemblyName, string typeName)
            {
                return Type.GetType(String.Format("{0}, {1}", typeName, assemblyName));
            }
        }

        /// <summary>
        ///     Convenience method for deserializing an object
        ///     BinaryFormatter allows you to serialize private member data unlike XmlSerializer
        /// </summary>
        public static T BinaryDeserialize<T>(string InputPath)
        {
            lock (_binarySerializationLock)
            {
                BinaryFormatter dat = new BinaryFormatter();
                dat.Binder = new MySerializationBinder(typeof(T));
                // allow shared read/write access in case file lock was not released even though serialization process has exited
                using (FileStream InputStream = new FileStream(InputPath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite))
                {
                    return (T)dat.Deserialize(InputStream);
                }
            }
        }

        /// <summary>
        ///     Convenience method for deserializing an object from XML
        /// </summary>
        public static T XMLDeserialize<T>(string inputPath)
        {
            XmlSerializer dat = new XmlSerializer(typeof(T));
            using (FileStream inputStream = new FileStream(inputPath, FileMode.Open, FileAccess.Read))
            {
                return (T)dat.Deserialize(inputStream);
            }
        }

        public static string GetMonoPath()
        {
            // We no longer package mono in the workflow directories
            // Assume it's already in the user's path
            return FindExePath("mono");
        }

        public static string FindExePath(string exe)
        {
            exe = Environment.ExpandEnvironmentVariables(exe);
            if (!File.Exists(exe))
            {
                if (Path.GetDirectoryName(exe) == String.Empty)
                {
                    foreach (string test in (Environment.GetEnvironmentVariable("PATH") ?? "").Split(':'))
                    {
                        string path = test.Trim();
                        if (!String.IsNullOrEmpty(path) && File.Exists(path = Path.Combine(path, exe)))
                            return Path.GetFullPath(path);
                    }
                }
                throw new FileNotFoundException(new FileNotFoundException().Message, exe);
            }
            return Path.GetFullPath(exe);
        }


        /// <summary>
        ///     Deletes only files that exist
        ///     return Exception if there were problems
        /// </summary>
        public static Exception SafeDelete(string filePath)
        {
            if (string.IsNullOrEmpty(filePath)) return null;
            if (File.Exists(filePath))
            {
                try
                {
                    File.Delete(filePath);
                }
                catch (Exception e)
                {
                    // Note: This shouldn't happen.  It could happen if (a) some other process has the file open and locked,
                    // or (b) if some other process deleted the file in between our call to File.Exists() and our call to 
                    // File.Delete().  Either way, we complain and then carry on.
                    Console.Error.WriteLine("Unable to delete file at {0}. Reason: {1}", filePath, e.ToString());
                    return e;
                }
            }
            return null;
        }

        /// <summary>
        ///     Deletes destination file before move
        /// </summary>
        public static void Move(string sourcePath, string destinationPath)
        {
            Directory.CreateDirectory(Path.GetDirectoryName(destinationPath));
            // note File.Delete does not throw an exception if the file does not exist
            File.Delete(destinationPath);
            File.Move(sourcePath, destinationPath);
        }

        /// <summary>
        ///     Deletes destination file before move
        /// </summary>
        public static void MoveFolder(string sourcePath, string destinationPath)
        {
            // note File.Delete does not throw an exception if the file does not exist
            if (Directory.Exists(destinationPath))
                Directory.Delete(destinationPath, true);
            Directory.CreateDirectory(Path.GetDirectoryName(destinationPath));
            Directory.Move(sourcePath, destinationPath);
        }

        public static List<Exception> SafeDelete(List<string> tempFilePaths)
        {
            var exceptions = new List<Exception>();
            foreach (string filePath in tempFilePaths)
            {
                var exception = SafeDelete(filePath);
                if (exception != null) exceptions.Add(exception);
            }
            return exceptions;
        }

        /// <summary>
        ///     Deletes only files that exist and only if IsasConfiguration does not have RetainTempFiles set
        ///     return Exception if there were problems
        /// </summary>
        public static Exception SafeDeleteTemp(string filePath)
        {
            if (IsasConfiguration.GetConfiguration().RetainTempFiles) return null;
            return SafeDelete(filePath);
        }
        public static List<Exception> SafeDeleteTemp(List<string> filePaths)
        {
            if (IsasConfiguration.GetConfiguration().RetainTempFiles) return null;
            return SafeDelete(filePaths);
        }

        /// <summary>
        /// Delete all files with a given suffix in one directory
        /// </summary>
        public static void SafeDeleteFilesBySuffix(string directory, string suffix)
        {
            foreach (FileInfo f in new DirectoryInfo(directory).GetFiles("*" + suffix))
            {
                Utilities.SafeDelete(f.FullName);
            }
        }

        /// <summary>
        /// Delete all files with a given suffix in one directory if retainTempFiles is not set
        /// </summary>
        public static void SafeDeleteTempFilesBySuffix(string directory, string suffix)
        {
            if (!IsasConfiguration.GetConfiguration().RetainTempFiles)
                SafeDeleteFilesBySuffix(directory, suffix);
        }

        public static void DeleteFolder(string path)
        {
            if (Directory.Exists(path))
                Directory.Delete(path, true);
        }

        /// <summary>
        ///     Delete a directory tree, catching and returning any exceptions.
        /// </summary>
        public static Exception SafeDeleteFolder(string filePath)
        {
            if (string.IsNullOrEmpty(filePath)) return null;
            if (!Directory.Exists(filePath)) return null;
            try
            {
                DirectoryInfo dirInfo = new DirectoryInfo(filePath);
                foreach (FileInfo file in dirInfo.GetFiles()) file.Delete();
                foreach (DirectoryInfo subDirectory in dirInfo.GetDirectories())
                    SafeDeleteFolder(subDirectory.FullName);
                Directory.Delete(filePath, true);
            }
            catch (Exception e)
            {
                // Note: This shouldn't happen.  It could happen if (a) some other process has the file open and locked,
                // or (b) if some other process deleted the file in between our call to File.Exists() and our call to 
                // File.Delete().  Either way, we complain and then carry on.
                Console.Error.WriteLine("Unable to delete file at {0}. Reason: {1}", filePath, e.ToString());
                return e;
            }
            return null;
        }

        /// <summary>
        ///     Delete a directory tree only if RetainTempFiles is not set in IsasConfiguration
        ///     catching and returning any exceptions.
        /// </summary>
        public static Exception SafeDeleteTempFolder(string filePath)
        {
            if (IsasConfiguration.GetConfiguration().RetainTempFiles) return null;
            return SafeDeleteFolder(filePath);
        }


        public static string GetReverseComplement(string dna)
        {
            StringBuilder result = new StringBuilder();
            for (int charIndex = dna.Length - 1; charIndex >= 0; charIndex--)
            {
                result.Append(GetComplement(dna[charIndex]));
            }
            return result.ToString();
        }

        /// <summary>
        /// Light-weight method to complement a single character
        /// </summary>
        /// <param name="nt"></param>
        /// <returns></returns>
        public static char GetComplement(char nt)
        {
            switch (nt)
            {
                case 'A':
                    return 'T';
                case 'C':
                    return 'G';
                case 'G':
                    return 'C';
                case 'T':
                    return 'A';
                case 'a':
                    return 't';
                case 'c':
                    return 'g';
                case 'g':
                    return 'c';
                case 't':
                    return 'a';
                default:
                    return nt;
            }
        }

        /// <summary>
        /// Return the platform-aware executable path for a given root filename. A simple sanity
        /// check is also performed to make sure that the file exists and if not an exception is thrown
        /// </summary>
        public static string GetExecutablePath(string fileName, Dictionary<string, string> customPaths, string workerRelativeBinDirectory = "")
        {
            string executableDirectory = Path.Combine(Utilities.GetAssemblyFolder(typeof(Utilities)), workerRelativeBinDirectory);
            return GetExecutablePath(fileName, executableDirectory, customPaths);
        }

        /// <summary>
        /// Return the platform-aware executable path for a given file name. A simple sanity
        /// check is also performed to make sure that the file exists and if not an exception is thrown
        /// </summary>
        public static string GetExecutablePath(string fileName, string executableDirectory, Dictionary<string, string> customPaths)
        {
            string lowerName = fileName.ToLowerInvariant();
            string lowerNameRoot = lowerName.ReplaceEnd(".exe", "");
            if (customPaths != null &&
                (customPaths.ContainsKey(lowerName) ||
                 customPaths.ContainsKey(lowerNameRoot)))
            {
                string overridePath = customPaths.ContainsKey(lowerName) ? customPaths[lowerName] : customPaths[lowerNameRoot];

                if (!File.Exists(overridePath))
                {
                    throw new Exception(string.Format("ERROR: Could not find the user-specified custom executable path location for the {0} executable at '{1}'.", fileName, overridePath));
                }
                Console.WriteLine("* Warning: Overriding executable path to '{0}'.  Custom executable paths are not officially supported and may lead to errors now or later in the workflow.",
                    overridePath);
                return overridePath;
            }
            string executablePath = Path.Combine(executableDirectory, fileName);

            if (!IsThisMono() && !executablePath.EndsWith(".exe") && !executablePath.EndsWith(".jar")) executablePath += ".exe";

            if (!File.Exists(executablePath))
            {
                throw new Exception(string.Format("ERROR: Could not find the {0} executable ({1}).", fileName, executablePath));
            }
            return executablePath;
        }

        public static string CheckForExecutableInPath(string executableName)
        {
            string[] paths = Environment.GetEnvironmentVariable("PATH").Split(Path.PathSeparator);
            foreach (string candidatePath in paths)
            {
                if (candidatePath.Any(pathChar => Path.GetInvalidPathChars().Contains(pathChar)))
                    continue;

                string fullCandidatePath = Path.Combine(candidatePath, executableName);
                if (File.Exists(fullCandidatePath))
                {
                    return fullCandidatePath;
                }
            }
            throw new Exception("Could not find executable " + executableName + " in path");
        }

        private static IWorkManager GetWorkManager(IsasConfiguration config, ErrorHandler log, ErrorHandler error)
        {
            ILogger logger = new ErrorLogger(log, error);
            LocalWorkManager local = new LocalWorkManager(logger, null, 0, config.MaximumMemoryGB, config.MaximumHoursPerProcess);
            IWorkManager manager = local;
            if (config.UseCluster)
            {
                manager = new SGEWorkManager(logger, local,
                    config.ClusterQueueName, config.ClusterParallelEnvironmentName);
            }
            return manager;
        }

        /// <summary>
        /// Run all tasks (on the cluster or locally), controlling parallelism based on the task requirements.
        /// </summary>
        /// <param name="tasks">All tasks to run</param>
        /// <param name="taskRequirements">Resource requirements for a single task</param>
        /// <param name="maxAttempts"> >1 means retry failed tasks (with less parallelism)</param>
        /// <param name="throwExceptionOnError">If set to False, continue with other tasks after one fails. Do not throw</param>
        /// <returns>True if all tasks finished successfully</returns>
        public static bool DoWorkParallel(IsasConfiguration config, ErrorHandler log, ErrorHandler error,
            List<UnitOfWork> tasks, TaskResourceRequirements taskRequirements, int maxAttempts = 1,
            bool throwExceptionOnError = true)
        {
            IWorkManager manager = GetWorkManager(config, log, error);
            return manager.DoWorkParallel(tasks, taskRequirements,
                maxAttempts: maxAttempts, throwExceptionOnError: throwExceptionOnError);
        }

        /// <summary>
        /// launches a single, single-threaded job locally
        /// </summary>
        public static void DoWorkSingleThread(IsasConfiguration config, ErrorHandler log, ErrorHandler error, UnitOfWork job,
            bool checkReturnCode = true, bool redirectOutput = true, bool throwExceptionOnError = true)
        {
            DoWorkParallelThreads(config, log, error, new List<UnitOfWork> { job }, 1, checkReturnCode: checkReturnCode,
                redirectOutput: redirectOutput, throwExceptionOnError: throwExceptionOnError);
        }

        /// <summary>
        ///     Handle work by launching up to MaximumThreads child threads locally.
        ///     Creates a estimated TaskResourceRequirements based on maximumThreads and available resources in this machine.
        ///     *** When possible set explicit TaskResourceRequirements *** 
        /// </summary>
        /// <returns>True if all tasks finished successfully</returns>
        public static bool DoWorkParallelThreads(IsasConfiguration config, ErrorHandler log, ErrorHandler error,
            List<UnitOfWork> tasks, int maximumThreads,
            bool checkReturnCode = true, bool redirectOutput = true, bool throwExceptionOnError = true)
        {
            // Estimated requirements object to attach to jobs, since none was supplied.
            // We still use the maximumThreads parameter do decide how many jobs to actually launch
            TaskResourceRequirements reqs = AvailableResources.GetMaxParallelThreadsRequirements(maximumThreads);
            return DoWorkParallelThreads(config, log, error, tasks, reqs, maximumThreads,
                    checkReturnCode, redirectOutput, throwExceptionOnError);
        }

        /// <summary>
        /// Handle work by launching threads locally, considering the given task requirements
        /// If maximumThreads >0, it controls the number of different jobs started at once (and each may receive more processors).
        /// </summary>
        /// <returns>True if all tasks finished successfully</returns>
        public static bool DoWorkParallelThreads(IsasConfiguration config, ErrorHandler log, ErrorHandler error,
            List<UnitOfWork> tasks, TaskResourceRequirements taskRequirements, int maximumThreads = 0,
            bool checkReturnCode = true, bool redirectOutput = true, bool throwExceptionOnError = true)
        {
            IWorkManager manager = GetWorkManager(config, log, error);
            return manager.DoWorkParallelLocal(tasks, taskRequirements,
                maximumThreads, checkReturnCode, redirectOutput, throwExceptionOnError);
        }

        /// <summary>
        /// A workaround for potential hang in Parallel.ForEach when using mono 4.0.2
        /// https://bugzilla.xamarin.com/show_bug.cgi?id=35297
        /// </summary>
        public static void DoWorkParallelThreads(List<ThreadStart> tasks, int? maxNumThreads = null)
        {
            if (maxNumThreads == null)
                maxNumThreads = MachineInfo.GetPhysicalProcessorCoreCount();

            // sanity check: do we have any jobs?
            if (tasks.Count == 0) return;

            // we originally used TPL, but we encountered unexpected behavior from MaxDegreeOfParallelism
            // With MaxDegreeOfParallelism, it seemed like if it was set to 10 and each job used 50% CPU,
            // it would launch around 20 threads. ThreadPool was better, but also tried to be smart about
            // launching threads. Using standard threads and a semaphore yielded the desired behavior.

            // run our tasks
            var threadPool = new Semaphore(maxNumThreads.GetValueOrDefault(), maxNumThreads.GetValueOrDefault());
            var doneEvent = new AutoResetEvent(false);
            var tasksRemaining = tasks.Count;
            foreach (ThreadStart task in tasks)
            {
                threadPool.WaitOne();

                Thread taskThread = new Thread(() =>
                {
                    task.Invoke();
                    if (Interlocked.Decrement(ref tasksRemaining) <= 0)
                        doneEvent.Set();
                    threadPool.Release();
                });
                taskThread.Start();
            }

            doneEvent.WaitOne();
        }

        // function to read text pos files as a last resort
        public static FloatPoint[] LoadPosFile(string locsPath)
        {
            List<FloatPoint> returnVals = new List<FloatPoint>();

            char[] splitchars = new[] { ' ', '\t' };
            using (StreamReader reader = new StreamReader(locsPath))
            {
                while (!reader.EndOfStream)
                {
                    string line = reader.ReadLine();
                    string[] bits = line.Split(splitchars, StringSplitOptions.RemoveEmptyEntries);
                    FloatPoint point = new FloatPoint { X = float.Parse(bits[0]), Y = float.Parse(bits[1]) };
                    returnVals.Add(point);
                }
            }

            return returnVals.ToArray();
        }

        public static void CreateSymbolicLinkOrCopy(string linkPath, string sourcePath)
        {
            if (!Path.IsPathRooted(sourcePath))
            {
                if (Directory.Exists(sourcePath))
                {
                    sourcePath = Path.Combine(linkPath, sourcePath);
                    new DirectoryLocation(sourcePath).CreateRelativeSymlinkOrCopy(new DirectoryLocation(linkPath));
                }
                else
                {
                    sourcePath = Path.Combine(Path.GetDirectoryName(linkPath), sourcePath);
                    new FileLocation(sourcePath).CreateRelativeSymlinkOrCopy(new FileLocation(linkPath));
                }
            }
            else if (Directory.Exists(sourcePath))
            {
                new DirectoryLocation(sourcePath).CreateAbsoluteSymlinkOrCopy(new DirectoryLocation(linkPath));
            }
            else
            {
                new FileLocation(sourcePath).CreateAbsoluteSymlinkAt(new FileLocation(linkPath));
            }
        }

        /// <summary>
        /// The real path to this file/directory after following any symlinks and removing any intermediate relative paths (e.g. "/../" or "/./")
        /// The actual file/directory does not need to exist
        /// NOTE: This method is deprecated. You should use the FullNameCanonical properties of IDirectoryLocation and IFileLocation
        /// </summary>
        public static string GetCanonicalPath(string path)
        {
            if (File.Exists(path))
            {
                return new FileLocation(path).FullNameCanonical;
            }
            else if (Directory.Exists(path))
            {
                return new DirectoryLocation(path).FullNameCanonical;
            }
            else
            {
                var parentDirectory = new DirectoryLocation(Path.GetDirectoryName(path));
                return Path.Combine(parentDirectory.FullNameCanonical, Path.GetFileName(path));
            }
        }

        /// <summary>
        ///     Waits until a lock has been released on the specified file.
        /// </summary>
        public static void WaitForFileUnlock(string filename, ErrorHandler error)
        {
            DateTime startTime = DateTime.Now;

            if (!File.Exists(filename)) return;
            while (true)
            {
                bool isLocked = false;

                // test if the file is currently locked
                try
                {
                    using (new FileStream(filename, FileMode.Open, FileAccess.Read, FileShare.None))
                    {
                        // deliberately do nothing
                    }
                }
                catch (IOException)
                {
                    isLocked = true;
                }

                // file is no longer locked
                if (!isLocked) break;

                if ((DateTime.Now - startTime).Seconds > 300)
                {
                    // We're done waiting; throw an exception!
                    throw new Exception(string.Format("Error: Gave up on waiting for unlock of file {0}", filename));
                }

                // perhaps we need to sleep 5 seconds at a time
                if (error != null) error(string.Format("Waiting 5 seconds for a file unlock on {0}", filename));
                Thread.Sleep(5000);
            }
        }

        public static string GetAssemblyFolder(Type assemblyType)
        {
            Assembly assembly = Assembly.GetAssembly(assemblyType);

            try
            {
                if (assembly != null)
                {
                    if (!string.IsNullOrEmpty(assembly.CodeBase))
                    {
                        string assemblyPath = new Uri(assembly.CodeBase).LocalPath;
                        return new FileInfo(assemblyPath).Directory.FullName;
                    }
                    Path.GetDirectoryName(assembly.Location);
                }
            }
            catch // (Exception ex)
            {
            }

            return string.Empty;
        }

        public static void UncompressFile(string sourcePath, string targetPath)
        {
            using (ZInputStream zstr = new ZInputStream(sourcePath))
            using (BinaryReader reader = new BinaryReader(zstr))
            using (FileStream outputStream = new FileStream(targetPath, FileMode.Create, FileAccess.Write))
            using (BinaryWriter writer = new BinaryWriter(outputStream))
            {
                while (true)
                {
                    byte[] buffer = reader.ReadBytes(1024 * 1024);
                    if (buffer == null || buffer.Length == 0) break;
                    writer.Write(buffer);
                }
            }
        }

        /// <summary>
        ///     Takes an input command line and merges in arbitrary options. Returns the updated command line
        ///     - If the option already exists in the command line, then the value is overridden.
        ///     - If the option does not exist in the command line, then the option and value are added after any other updated options
        ///     - If there are no other updated options then we insert at end or beginning which is controlled by insertAtEnd parameter
        ///     - Additional option #foo means: remove option foo if it's in the command-line now, otherwise do nothing.
        ///       (Remove "-foo" or "--foo" or "-foo bar" or "--foo bar")
        ///     "-param -4","--param -4" "--param=-4" and "-p-4" are all supported. 
        /// </summary>
        public static string MergeCommandLineOptions(string commandLineWithOptions, string moreOptions, bool insertAtEnd = false)
        {
            if (string.IsNullOrEmpty(moreOptions))
                return commandLineWithOptions;

            //parse the existing command line and the new options to generate list of option key/value pairs for each
            string beforeFirstOption;
            List<KeyValuePair<string, string>> updatedOptions = GetCommandOptions(commandLineWithOptions, out beforeFirstOption);
            string beforeNewOptions;
            List<KeyValuePair<string, string>> newOptions = GetCommandOptions(moreOptions, out beforeNewOptions);

            //validate the new options
            if (!string.IsNullOrEmpty(beforeNewOptions.Trim()))
                throw new ArgumentException(string.Format("Unknown options format {0}", moreOptions));

            //special case, the last option value may include additional arguments to the command that are not part of that option
            //if we are updating this last option we can try to infer the number of white space delimited values the option has (usually 0 or 1) from the new values
            string afterLastOption = "";
            foreach (KeyValuePair<string, string> option in newOptions)
            {
                if (option.Key == updatedOptions.Last().Key)
                {
                    string values = GetValues(updatedOptions.Last().Value, option.Value, out afterLastOption);
                    updatedOptions[updatedOptions.Count - 1] = new KeyValuePair<string, string>(option.Key, values);
                    break;
                }
            }

            //another special case, if we find a option to update make sure we do that first so we have a point of reference for inserting new options
            //inserting new options after this reference point has a better chance of being correct than arbitrarily choosing the beginning or the end of the command for insertion
            List<KeyValuePair<string, string>> prioritizedNewOptions = new List<KeyValuePair<string, string>>(newOptions);
            foreach (KeyValuePair<string, string> option in newOptions)
            {
                if (updatedOptions.Exists(kvp => option.Key == kvp.Key))
                {
                    prioritizedNewOptions.Remove(option);
                    prioritizedNewOptions.Insert(0, option);
                    break;
                }
            }

            //now we finally do the updating
            UpdateCommandOptions(updatedOptions, prioritizedNewOptions, insertAtEnd);

            //build up the final command that includes the update options
            StringBuilder updatedCommand = new StringBuilder(beforeFirstOption);
            foreach (KeyValuePair<string, string> option in updatedOptions)
            {
                updatedCommand.AppendWithSpace(option.Key);
                updatedCommand.Append(option.Value);
            }
            updatedCommand.AppendWithSpace(afterLastOption);
            return updatedCommand.ToString();
        }

        public static void UpdateCommandOptions(List<KeyValuePair<string, string>> originalOptions, List<KeyValuePair<string, string>> newOptions, bool insertAtEnd)
        {
            int lastUpdatedOptionIndex = -1;
            foreach (KeyValuePair<string, string> newOption in newOptions)
            {
                // Handle option removals:
                if (newOption.Key.StartsWith("#"))
                {
                    List<KeyValuePair<string, string>> removals = new List<KeyValuePair<string, string>>();
                    for (int index = 0; index < originalOptions.Count; index++)
                    {
                        var option = originalOptions[index];
                        if (option.Key.TrimStart('-') == newOption.Key.Substring(1))
                        {
                            removals.Add(option);
                            if (lastUpdatedOptionIndex >= index) lastUpdatedOptionIndex--;
                        }
                    }
                    foreach (var removal in removals) originalOptions.Remove(removal);
                    continue;
                }

                // Handle option overrides and insertions:
                int newOptionIndex = originalOptions.FindIndex(kvp => kvp.Key == newOption.Key);
                if (newOptionIndex != -1)
                {
                    originalOptions[newOptionIndex] = newOption;
                }
                else if (lastUpdatedOptionIndex != -1)
                    originalOptions.Insert(lastUpdatedOptionIndex + 1, newOption);
                else if (insertAtEnd)
                    originalOptions.Add(newOption);
                else
                    originalOptions.Insert(0, newOption);

                lastUpdatedOptionIndex = originalOptions.IndexOf(newOption);
            }
        }

        //count the number of values in the expectedValues string and grab the corresponding number of values from the currentValues string.
        //anything left over goes in the remainingValues string
        public static string GetValues(string currentValues, string expectedValues, out string remainingValues)
        {
            //a value can start at the beginning of the string with an optional = character and space --> \A=?\s*)
            //or simply by having space --> \s
            //the actual value is the a sequence of one or more non-whitespace characters --> \S+
            //finally the value must be followed by more space or the end of the string --> \s|\z
            Regex optionRegex = new Regex(@"((?:\A=?\s*|\s)\s*)(\S+)(?=\s|\z)");
            int expectedValuesCount = optionRegex.Matches(expectedValues).Count;
            MatchCollection matches = optionRegex.Matches(currentValues);
            Match lastValue = matches[expectedValuesCount - 1];
            int indexRemainingValues = lastValue.Index + lastValue.Length;
            string values = currentValues.Substring(0, indexRemainingValues);
            remainingValues = currentValues.Substring(indexRemainingValues);
            return values;
        }

        //split the command into option key/value pairs and separately save anything that comes before the first option 
        public static List<KeyValuePair<string, string>> GetCommandOptions(string command, out string beforeFirstOption)
        {
            //this regex is the meat of this feature. We want to capture all valid options from the command
            //each option must happen after white space or at the beginning of the command
            //the option must start with a hyphen and then may contain more hyphens only if they are followed by a non-digit
            //the end of the option is identified by an = character or whitespace or a hyphen followed by a digit or the end of the string 
            Regex optionRegex = new Regex(@"(?<=\A|\s)(((--?)|(#))[a-zA-Z](?:\-(?=\D)|[_a-zA-Z])*)(?=[=\s]|\-?\d|\z)");
            List<KeyValuePair<string, string>> options = new List<KeyValuePair<string, string>>();
            Match m = optionRegex.Match(command);
            if (m.Success)
                beforeFirstOption = command.Substring(0, m.Index);
            else
                beforeFirstOption = command;
            while (m.Success)
            {
                int valueStartIndex = m.Index + m.Length;
                Match next = m.NextMatch();
                string value;
                if (next.Success)
                    value = command.Substring(valueStartIndex, next.Index - valueStartIndex);
                else
                    value = command.Substring(valueStartIndex);
                options.Add(new KeyValuePair<string, string>(m.Value, value));
                m = m.NextMatch();
            }
            return options;
        }

        //add white space if necessary
        private static void AppendWithSpace(this StringBuilder start, string end)
        {
            if (start.Length > 0 && !Char.IsWhiteSpace(start[start.Length - 1]) && !string.IsNullOrEmpty(end) && !Char.IsWhiteSpace(end[0]))
                start.Append(" ");
            start.Append(end);
        }

        public static void DirectoryCopy(string sourceDirName, string destDirName, bool copySubDirs, bool overWrite = true)
        {
            // Get the subdirectories for the specified directory.
            DirectoryInfo dir = new DirectoryInfo(sourceDirName);
            DirectoryInfo[] dirs = dir.GetDirectories();

            if (!dir.Exists)
            {
                throw new DirectoryNotFoundException(
                    "Source directory does not exist or could not be found: "
                    + sourceDirName);
            }

            // If the destination directory doesn't exist, create it. 
            if (!Directory.Exists(destDirName))
            {
                Directory.CreateDirectory(destDirName);
            }

            // Get the files in the directory and copy them to the new location.
            FileInfo[] files = dir.GetFiles();
            foreach (FileInfo file in files)
            {
                string temppath = Path.Combine(destDirName, file.Name);
                file.CopyTo(temppath, overWrite);
            }

            // If copying subdirectories, copy them and their contents to new location. 
            if (copySubDirs)
            {
                foreach (DirectoryInfo subdir in dirs)
                {
                    string temppath = Path.Combine(destDirName, subdir.Name);
                    DirectoryCopy(subdir.FullName, temppath, copySubDirs, overWrite: overWrite);
                }
            }
        }

        public static string GetRelativePath(string filespec, string folder)
        {
            Uri pathUri = new Uri(filespec);
            // Folders must end in a slash
            if (!folder.EndsWith(Path.DirectorySeparatorChar.ToString()))
            {
                folder += Path.DirectorySeparatorChar;
            }
            Uri folderUri = new Uri(folder);
            return Uri.UnescapeDataString(folderUri.MakeRelativeUri(pathUri).ToString().Replace('/', Path.DirectorySeparatorChar));
        }
    }

    public class TaskResourceRequirements
    {
        /// Static parameters descripting available resources. Set once per Isas run.
        /// todo refactor this to separate per process requirements class and global singleton scheduler object.
        // Single node execution
        private static int NumCores { get { return AvailableResources.NumCores; } }
        private static double AvailableMemory { get { return AvailableResources.AvailableMemory; } }
        // SGE execution
        public static int CoresPerSGESlot { get { return AvailableResources.CoresPerSGESlot; } }
        public static float GigabytesPerSGESlot { get { return AvailableResources.GigabytesPerSGESlot; } }

        // Process requirements
        public readonly int MaxCores;
        public readonly double MemoryGBPerExtraCore;
        public readonly int MinCores;
        public readonly double MinMemoryGB;

        /// <summary>
        /// Set static machine resource limits, controling future scheduled jobs
        /// </summary>
        public static void InitializeAvailableResources(int MaximumThreadCount, IsasConfiguration config)
        {
            AvailableResources.InitializeAvailableResources(MaximumThreadCount, config.MaximumMemoryGB,
                config.ClusterCoresPerSlot, config.ClusterGigabytesPerSlot);
        }

        /// <summary>
        /// Requirements for single-threaded job with up to 4GB of RAM
        /// </summary>
        public static TaskResourceRequirements GetDefaultSingleThreadedTaskResourceRequirements()
        {
            return new TaskResourceRequirements(1, 4.0);
        }

        /// <summary>
        /// Requirements for a single task.
        /// Used to control degree of parallelism on a single node or SGE parameters.
        /// Used for jobs that allow different number of cores to be used (processes / threads).
        /// </summary>
        /// <param name="minCores">Minimum number of cores needed for one job</param>
        /// <param name="maxCores">Maximum number of cores useful for one job</param>
        /// <param name="minMemoryGB">RAM needed when run with minCores</param>
        /// <param name="memoryGBPerExtraCore">RAM needed for every extra core above minCores</param>
        public TaskResourceRequirements(int minCores, int maxCores, double minMemoryGB, double memoryGBPerExtraCore)
        {
            this.MinCores = minCores;
            this.MaxCores = maxCores;
            this.MinMemoryGB = minMemoryGB;
            this.MemoryGBPerExtraCore = memoryGBPerExtraCore;
        }

        /// <param name="Cores">Number of cores used by this job (threads / processes)</param>
        /// <param name="memoryGB">Max RAM used by job</param>
        public TaskResourceRequirements(int Cores, double memoryGB)
        {
            this.MinCores = Cores;
            this.MaxCores = Cores;
            this.MinMemoryGB = memoryGB;
            this.MemoryGBPerExtraCore = 0.0;
        }

        private void GetParallelParameters(int numberTasks, out int numberCoresPerTask, out int maximumParallelTasks)
        {
            numberCoresPerTask = MinCores; //Always run at least 1 Job with its minCores.
            maximumParallelTasks = 1;
            for (int assignedCores = MinCores; assignedCores <= MaxCores; assignedCores += 1)
            {
                // How many jobs could run in parallel with this number cores?
                int acceptableParallelJobsByCores = NumCores / assignedCores;
                int acceptableParallelJobsByMemory =
                    (int)(AvailableMemory / (MinMemoryGB + (assignedCores - MinCores) * MemoryGBPerExtraCore));

                int acceptableParallelJobs = Math.Min(acceptableParallelJobsByCores,
                                                      acceptableParallelJobsByMemory);

                // Accept this numer of cores if we can run as many jobs as there are tasks or
                // as many as we could run with fewer cores.
                if (acceptableParallelJobs >= maximumParallelTasks)
                {
                    numberCoresPerTask = assignedCores;
                    maximumParallelTasks = Math.Min(acceptableParallelJobs, numberTasks);
                }
            }
        }

        /// <summary>
        /// Max. number of SGE slots this job could use
        /// </summary>
        public int GetMaxSlots()
        {
            int maxSlots = 1;

            // How many slots would it take to satisfy my wildest dreams
            // for cores?
            //
            maxSlots = Math.Max((int)Math.Ceiling(((float)MaxCores) / CoresPerSGESlot), maxSlots);

            // How many cores would it take to satisfy the maximal
            // memory request?
            //
            maxSlots = Math.Max((int)Math.Ceiling((MinMemoryGB + (MemoryGBPerExtraCore) * (MaxCores - MinCores)) / GigabytesPerSGESlot), maxSlots);



            // Put some sane limits on the request
            //
            return Math.Min(maxSlots, 100);
        }

        /// <summary>
        /// Min. number of SGE slots this job needs
        /// </summary>
        public int GetMinSlots()
        {
            // Make sure we're requesting enough slots to cover both
            // our minimum cores and memory demands
            // 
            int minSlots = (int)Math.Ceiling(((float)MinCores / CoresPerSGESlot));
            minSlots = Math.Max((int)Math.Ceiling(MinMemoryGB / GigabytesPerSGESlot), minSlots);
            return minSlots;
        }


        /// <summary>
        /// When running on a single-node, how many cores to give to each job in this set
        /// </summary>
        /// <param name="numberTasks">Number of jobs in this set</param>
        public int GetNumberCoresPerTask(int numberTasks)
        {
            int cores;
            int parallelTasks;
            GetParallelParameters(numberTasks, out cores, out parallelTasks);
            return Math.Max(1, cores);
        }

        /// <summary>
        /// When running on a single-node, how many individual jobs to run at once
        /// </summary>
        /// <param name="numberTasks">Number of jobs in this set</param>
        public int GetMaximumParallelTasks(int numberTasks)
        {
            int cores;
            int parallelTasks;
            GetParallelParameters(numberTasks, out cores, out parallelTasks);
            return Math.Max(1, parallelTasks);
        }
    }

    /// <summary>
    /// This implementation of ILogger is to help the original, static DoWork methods use the 
    /// new IWorkManager interface.
    /// </summary>
    public class ErrorLogger : ILogger
    {
        private ErrorHandler error;
        private ErrorHandler log;

        public ErrorLogger(ErrorHandler log, ErrorHandler error)
        {
            this.error = error;
            this.log = log;
        }

        public void Dispose()
        {
            throw new NotImplementedException();
        }

        public void Log(LogEntry entry)
        {
            if (entry == null) return;

            if (entry.Severity >= LogEventType.Warning)
            {
                error(entry.Message);
            }
            else
            {
                log(entry.Message);
            }
        }
    }
}
