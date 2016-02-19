using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Illumina.SecondaryAnalysis;
using ILMNcommon.Common;

namespace Isas.Shared.Checkpointing
{
    public class CheckpointRunnerAsync : ICheckpointRunnerAsync
    {
        private readonly IDirectoryLocation _tempRepository;
        private readonly CheckpointManagerAsync _manager;
        private readonly ICheckpointSerializerAsync _serializer;
        private readonly List<string> _childStartCheckpoints;
        private readonly List<string> _childStopCheckpoints;
        private readonly ILogger _logger;
        private readonly bool _retainTempFiles;
        private readonly bool _forceSynchronous;
        private readonly List<Checkpoint> _parentCheckpoints;
        private bool _disposed;
        private Task _stopCheckpointTask;

        private CheckpointRunnerAsync(ILogger logger, IDirectoryLocation tempRepository, CheckpointManagerAsync manager, ICheckpointSerializerAsync serializer,
            bool retainTempFiles = false, bool forceSynchronous = false, List<string> childStartCheckpoints = null, List<string> childStopCheckpoints = null, List<Checkpoint> parentCheckpoints = null)
        {
            _tempRepository = tempRepository;
            _manager = manager;
            _serializer = serializer;
            _childStartCheckpoints = childStartCheckpoints ?? new List<string>();
            _childStopCheckpoints = childStopCheckpoints ?? new List<string>();
            _logger = logger;
            _retainTempFiles = retainTempFiles;
            _forceSynchronous = forceSynchronous;
            _parentCheckpoints = parentCheckpoints ?? new List<Checkpoint>();
        }

        private CheckpointHelper<T> CreateCheckpoint<T>(Checkpoint checkpoint, Func<T> func)
        {
            return new CheckpointHelper<T>(checkpoint, func, _logger, this);
        }

        public Task<TOut> RunCheckpointAsync<TOut>(string checkpointName, Func<TOut> func)
        {
            var checkpoint = _manager.CreateCheckpoint(checkpointName);
            return RunCheckpointAsync(checkpoint, func);
        }

        private Task<TOut> RunCheckpointAsync<TOut>(Checkpoint checkpoint, Func<TOut> func)
        {
            var checkpointHelper = CreateCheckpoint(checkpoint, func);
            return checkpointHelper
                .WrapWithBenchmark()
                .WrapWithSave()
                .RunOrLoad()
                .RunAsync();
        }

        private string FullCheckpoint(Checkpoint currentCheckpoint)
        {
            return string.Join(" --> ", GetParentCheckpoints(currentCheckpoint));
        }

        private string PrettyFullCheckpoint(Checkpoint currentCheckpoint)
        {
            if (_parentCheckpoints.Any())
                return $"nested checkpoint {FullCheckpoint(currentCheckpoint)}";
            return $"checkpoint {FullCheckpoint(currentCheckpoint)}";
        }

        public Task<TOutput> RunCheckpointAsync<TOutput>(
            string checkpointName,
            Func<ICheckpointRunnerAsync, IDirectoryLocation, TOutput> run)
        {
            var checkpoint = _manager.CreateCheckpoint(checkpointName);
            return RunCheckpointAsync(checkpoint, SupplyNestedCheckpoint(checkpoint, run));
        }

        private Func<TOutput> SupplyNestedCheckpoint<TOutput>(Checkpoint checkpoint, Func<ICheckpointRunnerAsync, IDirectoryLocation, TOutput> run)
        {
            return () =>
            {
                string startCheckpoint = null;
                var childStartCheckpoints = new List<string>();
                if (checkpoint.IsStartCheckpoint)
                {
                    startCheckpoint = _childStartCheckpoints.FirstOrDefault();
                    childStartCheckpoints = _childStartCheckpoints.Skip(1).ToList();
                }
                string stopCheckpoint = null;
                var childStopCheckpoints = new List<string>();
                if (checkpoint.IsStopCheckpoint)
                {
                    stopCheckpoint = _childStopCheckpoints.FirstOrDefault();
                    childStopCheckpoints = _childStopCheckpoints.Skip(1).ToList();
                }
                IDirectoryLocation tempDir = GetTempDirectory(checkpoint);
                var childSerializer = _serializer.GetChildSerializer(CheckpointManagerFactory.GetCheckpointFolder(tempDir));

                var childManager = _manager.CreateChildCheckpointerManager(startCheckpoint, stopCheckpoint);
                using (var childCheckpointer = GetChildCheckpointer(tempDir, childManager, childSerializer, childStartCheckpoints, childStopCheckpoints, GetParentCheckpoints(checkpoint)))
                {
                    if (!childStartCheckpoints.Any())
                        tempDir.CreateClean();
                    var results = run.Invoke(childCheckpointer, tempDir);
                    if (!_retainTempFiles)
                        tempDir.Delete();
                    return results;
                }
            };
        }

        private List<Checkpoint> GetParentCheckpoints(Checkpoint currentCheckpoint)
        {
            var parentCheckpoints = new List<Checkpoint>();
            parentCheckpoints.AddRange(_parentCheckpoints);
            parentCheckpoints.Add(currentCheckpoint);
            return parentCheckpoints;
        }

        private CheckpointRunnerAsync GetChildCheckpointer(IDirectoryLocation tempDir, CheckpointManagerAsync childManager, ICheckpointSerializerAsync childSerializer,
            List<string> childStartCheckpoints, List<string> childStopCheckpoints, List<Checkpoint> parentCheckpoints)
        {
            return new CheckpointRunnerAsync(_logger, tempDir, childManager, childSerializer, _retainTempFiles, _forceSynchronous, childStartCheckpoints, childStopCheckpoints, parentCheckpoints);
        }

        public Task<TOut> RunCheckpointAsync<TIn, TOut>(string key, Func<ICheckpointRunnerAsync, TIn, IDirectoryLocation, TOut> run,
    TIn input, ILoadingConvention<TIn, TOut> loadingConvention)
        {
            var checkpoint = _manager.CreateCheckpoint(key);
            var checkpointHelper = CreateCheckpoint(checkpoint,
                SupplyNestedCheckpoint(checkpoint, WrapWithLoadingConvention(checkpoint, run, input, loadingConvention)));
            return checkpointHelper.RunAsync();
        }

        private Func<ICheckpointRunnerAsync, IDirectoryLocation, TOutput> WrapWithLoadingConvention<TInput, TOutput>(Checkpoint checkpoint, Func<ICheckpointRunnerAsync, TInput, IDirectoryLocation, TOutput> run,
            TInput input,
            ILoadingConvention<TInput, TOutput> convention)
        {
            return (runner, tempDir) =>
            {
                bool hasBeenRunBefore = _serializer.CanLoad(checkpoint);
                if (!checkpoint.ShouldRun(hasBeenRunBefore) && !hasBeenRunBefore)
                {
                    return CreateCheckpoint(checkpoint, () =>
                    {
                        _logger.Warn($"Skipping {PrettyFullCheckpoint(checkpoint)} by loading results from the filesystem. " +
                                     $"Any results from this checkpoint that are generated only in memory will be unavailable to downstream steps and may cause unexpected behavior");
                        return convention.Load(input);
                    }).WrapWithSave().Func.Invoke();
                }
                return CreateCheckpoint(checkpoint, () =>
                {
                    var results = run.Invoke(runner, input, tempDir);
                    convention.Move(results, (source, destination) => source.MoveAndLink(destination));
                    return convention.Load(input);
                }).WrapWithBenchmark().WrapWithSave().RunOrLoad().Func.Invoke();
            };
        }

        public Task<TOutput> RunCheckpointAsync<TInput, TOutput>(
            string checkpointName,
            Func<TInput, IDirectoryLocation, TOutput> run,
            TInput input,
            ILoadingConvention<TInput, TOutput> convention)
        {
            return RunCheckpointAsync(checkpointName, (checkpointer, runInput, tempDir) => run(runInput, tempDir), input, convention);
        }

        public Task<TOutput> RunCheckpointAsync<TInput, TOutput>(
            string checkpointName,
            Func<TInput, IDirectoryLocation, TOutput> run,
            TInput input,
            INamingConvention<TOutput> convention)
        {
            var checkpoint = _manager.CreateCheckpoint(checkpointName);
            return RunCheckpointAsync(checkpoint, SupplyNestedCheckpoint(checkpoint, (nestedRunner, tempDir) =>
            {
                var results = run.Invoke(input, tempDir);
                return convention.Move(results, (source, destination) => source.MoveAndLink(destination));
            }));
        }

        private IDirectoryLocation GetTempDirectory(Checkpoint checkpoint)
        {
            return _tempRepository.GetDirectoryLocation($"Temp{checkpoint.Number:00}-{checkpoint.Name.RemoveWhiteSpace()}");
        }

        public void Dispose()
        {
            if (_disposed) return;
            _manager.Dispose();
            _disposed = true;
        }

        public TOut RunCheckpoint<TOut>(string key, Func<TOut> function)
        {
            return RunCheckpointAsync(key, function).GetAwaiter().GetResult();
        }

        public TOut RunCheckpoint<TOut>(string key, Func<IDirectoryLocation, TOut> function)
        {
            return RunCheckpointAsync(key, (checkpointer, tempDir) => function(tempDir)).GetAwaiter().GetResult();
        }

        public TOutput RunCheckpoint<TOutput>(string checkpointName, Func<ICheckpointRunner, IDirectoryLocation, TOutput> run)
        {
            return RunCheckpointAsync(checkpointName, run).GetAwaiter().GetResult();
        }

        public TOut RunCheckpoint<TIn, TOut>(string key, Func<TIn, IDirectoryLocation, TOut> run, TIn input, ILoadingConvention<TIn, TOut> loadingConvention)
        {
            return RunCheckpointAsync(key, run, input, loadingConvention).GetAwaiter().GetResult();
        }

        public TOut RunCheckpoint<TIn, TOut>(string key, Func<ICheckpointRunner, TIn, IDirectoryLocation, TOut> run, TIn input, ILoadingConvention<TIn, TOut> loadingConvention)
        {
            return RunCheckpointAsync(key, run, input, loadingConvention).GetAwaiter().GetResult();
        }

        public TOut RunCheckpoint<TIn, TOut>(string key, Func<TIn, IDirectoryLocation, TOut> run, TIn input, INamingConvention<TOut> namingConvention)
        {
            return RunCheckpointAsync(key, run, input, namingConvention).GetAwaiter().GetResult();
        }

        public static ICheckpointRunnerAsync Create(ICheckpointSerializerAsync serializer, ILogger logger, IDirectoryLocation tempDir,
            string startCheckpoint = null, string stopCheckpoint = null, bool retainTempFiles = false, bool forceSynchronous = false)
        {
            var startCheckpoints = string.IsNullOrEmpty(startCheckpoint) ? new List<string>() : startCheckpoint.Split('.').ToList();
            if (startCheckpoints.Any(string.IsNullOrWhiteSpace))
                throw new ArgumentException("Starting checkpoints cannot be empty");
            var invalidNumberedCheckpoints = startCheckpoints.Where(checkpointString =>
            {
                int checkpointNumber;
                if (int.TryParse(checkpointString, out checkpointNumber) && checkpointNumber <= 0)
                    return true;
                return false;
            });
            if (invalidNumberedCheckpoints.Any())
                throw new ArgumentException("Numbered starting checkpoints must be greater than or equal to 1");


            var stopCheckpoints = string.IsNullOrEmpty(stopCheckpoint) ? new List<string>() : stopCheckpoint.Split('.').ToList();
            if (stopCheckpoints.Any(string.IsNullOrWhiteSpace))
                throw new ArgumentException("Stopping checkpoints cannot be empty");
            invalidNumberedCheckpoints = stopCheckpoints.Where(checkpointString =>
            {
                int checkpointNumber;
                if (int.TryParse(checkpointString, out checkpointNumber) && checkpointNumber <= 0)
                    return true;
                return false;
            });
            if (invalidNumberedCheckpoints.Any())
                throw new ArgumentException("Numbered stopping checkpoints must be greater than or equal to 1");

            var manager = new CheckpointManagerAsync(logger, startCheckpoints.FirstOrDefault(), stopCheckpoints.FirstOrDefault());
            return new CheckpointRunnerAsync(logger, tempDir, manager, serializer, retainTempFiles, forceSynchronous, startCheckpoints.Skip(1).ToList(),
                stopCheckpoints.Skip(1).ToList());
        }

        private class CheckpointHelper<T>
        {
            private Checkpoint Checkpoint { get; }
            public Func<T> Func { get; }
            private readonly ILogger _logger;
            private readonly CheckpointRunnerAsync _runner;

            public CheckpointHelper(Checkpoint checkpoint, Func<T> func, ILogger logger, CheckpointRunnerAsync runner)
            {
                Checkpoint = checkpoint;
                Func = func;
                _logger = logger;
                _runner = runner;
            }

            public CheckpointHelper<T> WrapWithBenchmark()
            {
                return new CheckpointHelper<T>(Checkpoint, () =>
                {
                    _logger.Info($"Running {_runner.PrettyFullCheckpoint(Checkpoint)}");
                    Benchmark benchmark = new Benchmark();
                    var result = Func.Invoke();
                    _logger.Info("Elapsed time (step/time(sec)/name)\t{0}\t{1:F1}\t{2}",
                        Checkpoint.Number, benchmark.GetElapsedTime(), Checkpoint.Name);
                    return result;
                }, _logger, _runner);
            }

            public CheckpointHelper<T> WrapWithSave()
            {
                return new CheckpointHelper<T>(Checkpoint, () =>
                {
                    var result = Func.Invoke();
                    _runner._serializer.Save(Checkpoint, result);
                    return result;
                }, _logger, _runner);
            }

            public CheckpointHelper<T> RunOrLoad()
            {
                bool hasBeenRunBefore = _runner._serializer.CanLoad(Checkpoint);
                var func = Func;
                if (!Checkpoint.ShouldRun(hasBeenRunBefore))
                    func = () =>
                    {
                        _logger.Info($"Loading previously saved results for {_runner.PrettyFullCheckpoint(Checkpoint)}");
                        return _runner._serializer.Load<T>(Checkpoint);
                    };
                return new CheckpointHelper<T>(Checkpoint, func, _logger, _runner);
            }

            private CheckpointHelper<T> WrapWithStopCheckpointCheck()
            {
                return new CheckpointHelper<T>(Checkpoint, () =>
                {
                    var result = Func.Invoke();
                    if (Checkpoint.IsStopCheckpoint)
                    {
                        var message = $"Stop checkpoint {_runner.FullCheckpoint(Checkpoint)} finished";
                        _logger.Info(message);
                        throw new StopCheckpointFoundException(message);
                    }
                    return result;
                }, _logger, _runner);
            }

            private Task<T> RunAsyncInternal()
            {
                // handle any existing stop checkpoint
                if (_runner._stopCheckpointTask != null)
                {
                    Func<Task<T>> func = async () =>
                    {
                        await _runner._stopCheckpointTask;
                        // the above await should always throw a StopCheckpointFoundException
                        throw new ApplicationException("Unexpected program execution. A StopCheckpointFoundException was expected");
                    };
                    return func.Invoke();
                }

                Task<T> task;
                if (_runner._forceSynchronous)
                    task = Task.FromResult(Func.Invoke());
                else
                    task = Task.Run(Func);

                if (Checkpoint.IsStopCheckpoint)
                    _runner._stopCheckpointTask = task;
                return task;
            }

            public Task<T> RunAsync()
            {
                return WrapWithStopCheckpointCheck().RunAsyncInternal();
            }
        }
    }
}