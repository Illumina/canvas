using Drmaa;
using Illumina.SecondaryAnalysis;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Isas.Shared
{
    public interface IWorkManager
    {
        IDirectoryLocation LoggingFolder { get; }
        int MaximumThreadCount { get; }

        /// <summary>
        /// Run all tasks (on the cluster or locally), controlling parallelism based on the task requirements.
        /// </summary>
        /// <param name="tasks">All tasks to run</param>
        /// <param name="taskRequirements">Resource requirements for a single task</param>
        /// <param name="maxAttempts"> >1 means retry failed tasks (with less parallelism)</param>
        /// <param name="throwExceptionOnError">If set to False, continue with other tasks after one fails. Do not throw</param>
        /// <returns>True if all tasks finished successfully</returns>
        bool DoWorkParallel(List<UnitOfWork> tasks, TaskResourceRequirements taskRequirements,
            int maxAttempts = 1, bool throwExceptionOnError = true);

        /// <summary>
        /// Handle work by launching threads locally, considering the given task requirements
        /// If maximumThreads >0, it controls the number of different jobs started at once (and each may receive more processors).
        /// </summary>
        /// <returns>True if all tasks finished successfully</returns>
        bool DoWorkParallelLocal(List<UnitOfWork> tasks, TaskResourceRequirements taskRequirements, int maximumThreads = 0,
            bool checkReturnCode = true, bool redirectOutput = true, bool throwExceptionOnError = true);

        /// <summary>
        /// Executes the <paramref name="function"/> on each input in parallel, 
        /// using up to the maximum number of threads specified by the IWorkManager.
        /// Does not respect RAM limits. 
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="inputs"></param>
        /// <param name="function"></param>
        void ParallelForEach<T>(IEnumerable<T> inputs, Action<T> function);
    }

    public static class IWorkManagerExtensions
    {
        /// <summary>
        ///     Handle work by launching up to MaximumThreads child threads locally.
        ///     Creates a estimated TaskResourceRequirements based on maximumThreads and available resources in this machine.
        ///     *** When possible set explicit TaskResourceRequirements *** 
        /// </summary>
        public static bool DoWorkParallelThreads(this IWorkManager manager, List<UnitOfWork> tasks, int maximumThreads,
            bool checkReturnCode = true, bool redirectOutput = true, bool throwExceptionOnError = true)
        {
            // Estimated requirements object to attach to jobs, since none was supplied.
            // We still use the maximumThreads parameter do decide how many jobs to actually launch
            TaskResourceRequirements reqs = AvailableResources.GetMaxParallelThreadsRequirements(maximumThreads);
            return manager.DoWorkParallelLocal(tasks, reqs, maximumThreads,
                    checkReturnCode, redirectOutput, throwExceptionOnError);
        }

        /// <summary>
        ///     Handle work by launching up to MaximumThreadCount threads.
        ///     MaximumThreadCount can be specified in the SampleSheet (default is Environment.ProcessorCount)
        ///     Creates a estimated TaskResourceRequirements based on  MaximumThreadCount
        ///      and available resources in this machine.
        ///     *** When possible set explicit TaskResourceRequirements *** 
        /// </summary>
        public static void DoWorkParallelThreads(this IWorkManager manager, List<UnitOfWork> tasks,
            bool checkReturnCode = true, bool redirectOutput = true, bool throwExceptionOnError = true)
        {
            // Estimated requirements object to attach to jobs, since none was supplied.
            // We still use the maximumThreads parameter do decide how many jobs to actually launch
            TaskResourceRequirements reqs = AvailableResources.GetMaxParallelThreadsRequirements(manager.MaximumThreadCount);
            manager.DoWorkParallelLocal(tasks, reqs,
                    checkReturnCode: checkReturnCode, redirectOutput: redirectOutput,
                    throwExceptionOnError: throwExceptionOnError);
        }

        public static void DoWorkSingleThread(this IWorkManager manager,
            UnitOfWork job, bool checkReturnCode = true, bool redirectOutput = true, bool throwExceptionOnError = true)
        {
            TaskResourceRequirements reqs = AvailableResources.GetMaxParallelThreadsRequirements(1);
            manager.DoWorkParallelLocal(new List<UnitOfWork> { job }, reqs, 1,
                    checkReturnCode, redirectOutput, throwExceptionOnError);
        }

        public static UnitOfWork CreateUnitOfWork(this IWorkManager manager, string commandLine = null,
            string executablePath = null, string loggingStub = null)
        {
            return new UnitOfWork()
            {
                CommandLine = commandLine,
                ExecutablePath = executablePath,
                LoggingFolder = manager.LoggingFolder.FullName,
                LoggingStub = loggingStub
            };
        }

        public static void ParallelForEach<T>(this IWorkManager manager,
            SampleSet<T> inputs, Action<SampleInfo, T> function)
        {
            manager.ParallelForEach(inputs, kvp => function(kvp.Key, kvp.Value));
        }

        public static void ParallelForEach<T1, T2>(this IWorkManager manager,
            SampleSet<T1> inputs, SampleSet<T2> inputs2, Action<SampleInfo, T1, T2> function)
        {
            manager.ParallelForEach(
                inputs.Join(inputs2, (i1, i2) => new Tuple<T1, T2>(i1, i2)),
                kvp => function(kvp.Key, kvp.Value.Item1, kvp.Value.Item2));
        }

        public static void ParallelForEach<T1, T2, T3>(this IWorkManager manager,
            SampleSet<T1> inputs, SampleSet<T2> inputs2, SampleSet<T3> inputs3, Action<SampleInfo, T1, T2, T3> function)
        {
            manager.ParallelForEach(
                inputs.Join(inputs2, inputs3, (i1, i2, i3) => new Tuple<T1, T2, T3>(i1, i2, i3)),
                kvp => function(kvp.Key, kvp.Value.Item1, kvp.Value.Item2, kvp.Value.Item3));
        }
    }

    public class SGEWorkManager : IWorkManager
    {
        public IDirectoryLocation LoggingFolder => _localManager.LoggingFolder;
        public int MaximumThreadCount => _localManager.MaximumThreadCount;

        private readonly string _clusterQueueName;
        private readonly string _clusterParallelEnvironmentName;

        private readonly ILogger _logger;
        private readonly LocalWorkManager _localManager;

        public SGEWorkManager(ILogger logger, LocalWorkManager localManager,
            string queueName, string parallelEnvironmentName)
        {
            _logger = logger;
            _localManager = localManager;

            _clusterQueueName = queueName;
            _clusterParallelEnvironmentName = parallelEnvironmentName;
        }

        private void DoWorkOnCluster(List<UnitOfWork> tasks, TaskResourceRequirements taskRequirements, bool throwExceptionOnError)
        {
            using (Interop dInterop = new Interop(_clusterQueueName, _clusterParallelEnvironmentName,
                _logger.Info, _logger.Error))
            {
                dInterop.RunJobs(tasks, taskRequirements, throwExceptionOnError);
            }
        }

        public bool DoWorkParallel(List<UnitOfWork> tasks, TaskResourceRequirements taskRequirements, int maxAttempts = 1, bool throwExceptionOnError = true)
        {
            maxAttempts = Math.Max(1, maxAttempts);
            bool throwException = false;
            for (int attempt = 1; attempt <= maxAttempts && tasks.Count > 0; attempt++)
            {
                if (attempt == maxAttempts) throwException = throwExceptionOnError;
                _logger.Info($"DoWorkParallelThreads, attempt {attempt} of {maxAttempts}");

                DoWorkOnCluster(tasks, taskRequirements, throwException);

                tasks.RemoveAll(t => t.ExitCode == 0);

                // Reduce Cores/Task and increase memory after failure for retry
                taskRequirements = new TaskResourceRequirements(taskRequirements.MinCores,
                                                                taskRequirements.MinMemoryGB * 2);
            }
            return (tasks.Count == 0); // All done successfully
        }

        public bool DoWorkParallelLocal(List<UnitOfWork> tasks, TaskResourceRequirements taskRequirements,
            int maximumThreads = 0, bool checkReturnCode = true,
            bool redirectOutput = true, bool throwExceptionOnError = true)
        {
            return _localManager.DoWorkParallelLocal(tasks, taskRequirements, maximumThreads,
                checkReturnCode, redirectOutput, throwExceptionOnError);
        }

        public void ParallelForEach<T>(IEnumerable<T> inputs, Action<T> function)
        {
            _localManager.ParallelForEach(inputs, function);
        }
    }

    public class LocalWorkManager : IWorkManager
    {
        private readonly ILogger _logger;
        public int MaximumThreadCount { get; private set; }
        public IDirectoryLocation LoggingFolder { get; private set; }
        private float MaximumHoursPerProcess;
        private float MaximumMemoryGB;

        public LocalWorkManager(ILogger logger, IDirectoryLocation loggingFolder,
            int maximumThreadCount, float maximumMemoryGB, float maximumHoursPerProcess)
        {
            _logger = logger;
            if (loggingFolder != null)
            {
                LoggingFolder = loggingFolder;
                if (!LoggingFolder.Exists)
                    LoggingFolder.Create();
            }
            MaximumThreadCount = maximumThreadCount;
            MaximumMemoryGB = maximumMemoryGB;
            MaximumHoursPerProcess = maximumHoursPerProcess;

            AvailableResources.InitializeAvailableResources(maximumThreadCount, maximumMemoryGB, 0, 0);
        }

        public bool DoWorkParallel(List<UnitOfWork> tasks, TaskResourceRequirements taskRequirements, int maxAttempts = 1, bool throwExceptionOnError = true)
        {
            maxAttempts = Math.Max(1, maxAttempts);
            bool throwException = false;
            for (int attempt = 1; attempt <= maxAttempts && tasks.Count > 0; attempt++)
            {
                if (attempt == maxAttempts) throwException = throwExceptionOnError;
                _logger.Info($"DoWorkParallelThreads, attempt {attempt} of {maxAttempts}");

                DoWorkParallelLocal(tasks, taskRequirements, 0, true, true, throwException);
                tasks.RemoveAll(t => t.ExitCode == 0);

                // Reduce Cores/Task and increase memory after failure for retry
                taskRequirements = new TaskResourceRequirements(taskRequirements.MinCores,
                                                                taskRequirements.MinMemoryGB * 2);
            }
            return (tasks.Count == 0); // All done successfully
        }

        public bool DoWorkParallelLocal(List<UnitOfWork> tasks, TaskResourceRequirements taskRequirements,
            int maximumThreads = 0, bool checkReturnCode = true, bool redirectOutput = true, bool throwExceptionOnError = true)
        {
            if (maximumThreads == 0)
            {
                maximumThreads = taskRequirements.GetMaximumParallelTasks(tasks.Count);
            }
            int NumberCores = taskRequirements.GetNumberCoresPerTask(maximumThreads);
            foreach (UnitOfWork task in tasks) 
            {
                task.CommandLine = task.CommandLine.Replace(CommandLineTokens.NumProcessors, string.Format("{0}", NumberCores));
            }
            // todo pass TaskRequirements to the jobManager and use it there for monitoring etc.
            JobManager jobManager = new JobManager(MaximumMemoryGB, MaximumHoursPerProcess, _logger.Info, _logger.Error);
            jobManager.ProcessJobs(tasks, maximumThreads, redirectOutput, checkReturnCode, throwExceptionOnError);
            if (!checkReturnCode) return true;
            return tasks.All(t => t.ExitCode == 0); // All done successfully
        }

        public void ParallelForEach<T>(IEnumerable<T> inputs, Action<T> function)
        {
            var parOps = new ParallelOptions()
            {
                MaxDegreeOfParallelism = AvailableResources.NumCores
            };
            Parallel.ForEach(inputs, parOps, function);
        }

    }

    public static class AvailableResources
    {
        // Static parameters descripting available resources. Set once per Isas run.
        // todo refactor this to separate per process requirements class and global singleton scheduler object.

        // Single node execution
        public static int NumCores { get; private set; }
        public static double AvailableMemory { get; private set; }
        public static int MaximumThreadCount { get; private set; }

        // SGE execution
        public static int CoresPerSGESlot { get; private set; }
        public static float GigabytesPerSGESlot { get; private set; }

        /// <summary>
        /// Set static machine resource limits, controling future scheduled jobs
        /// </summary>
        public static void InitializeAvailableResources(int MaximumThreadCount, float maximumMemoryGB,
            int clusterCoresPerSlot, float clusterGigabytesPerSlot)
        {
            CoresPerSGESlot = clusterCoresPerSlot;
            GigabytesPerSGESlot = clusterGigabytesPerSlot;
            AvailableMemory = maximumMemoryGB;
            AvailableResources.MaximumThreadCount = MaximumThreadCount;
            NumCores = MaximumThreadCount > 0 ?
                    Math.Min(Environment.ProcessorCount, MaximumThreadCount) :
                    Environment.ProcessorCount;
        }

        /// <summary>
        /// Creates a (dummy) requirements object dividing the resources available on this machine into maximumParallelTasks equal pieces.
        /// Not used for actual scheduling at the moment, only to attach estimated requirements to jobs started through DoWorkParallelThreads
        /// which do not have an explicit TaskResourceRequirement attached.
        /// </summary>
        public static TaskResourceRequirements GetMaxParallelThreadsRequirements(int maximumParallelTasks)
        {
            if (maximumParallelTasks <= 0)
                maximumParallelTasks = NumCores;
            int requiredCores = Math.Max(1, (int)Math.Floor(NumCores / (double)maximumParallelTasks));
            double requiredMemory = Math.Max(.5, AvailableMemory - .5) / maximumParallelTasks;
            return new TaskResourceRequirements(requiredCores, requiredMemory);
        }
    }
}
