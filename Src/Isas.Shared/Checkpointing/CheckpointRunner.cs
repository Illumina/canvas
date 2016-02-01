using System;
using System.Linq;
using Illumina.SecondaryAnalysis;
using ILMNcommon.Common;

namespace Isas.Shared.Checkpointing
{
    public class CheckpointRunner : ICheckpointRunner
    {
        private readonly ILogger _logger;
        private readonly IDirectoryLocation _placeForTempFolders;
        private readonly ICheckpointSerializer _serializer;
        private readonly bool _retainTemps;
        private readonly ICheckpointManager _manager;
        private bool _disposed;

        public CheckpointRunner(ILogger logger, IDirectoryLocation placeForTempFolders, ICheckpointManager manager, ICheckpointSerializer serializer, bool retainTemps = false)
        {
            _logger = logger;

            _placeForTempFolders = placeForTempFolders;
            if (!_placeForTempFolders.Exists)
                _placeForTempFolders.Create();

            _manager = manager;
            _serializer = serializer;
            _retainTemps = retainTemps;
        }

        public TOut RunCheckpoint<TOut>(string key, Func<TOut> function)
        {
            _manager.BeginCheckpoint(key);
            TOut output;
            bool hasExistingSerializedOutput = HasExistingSerializedOutput();
            CheckpointStatus status = _manager.GetCheckpointStatus(hasExistingSerializedOutput);
            if (status.ShouldRun())
            {
                string checkpoint = GetCheckpointName();
                string checkpointNumber = GetCheckpointNumber();
                _logger.Info("Running checkpoint {0}: {1}", checkpointNumber, checkpoint);
                Benchmark benchmark = new Benchmark();
                output = function();
                _logger.Info("Elapsed time (step/time(sec)/name)\t{0}\t{1:F1}\t{2}",
                    checkpointNumber, benchmark.GetElapsedTime(), checkpoint);
                Save(output);
            }
            else
            {
                output = Load<TOut>();
            }
            _manager.EndCheckpoint();
            return output;
        }

        public TOutput RunCheckpoint<TInput, TOutput>(
            string key,
            Func<TInput, IDirectoryLocation, TOutput> run,
            TInput input,
            INamingConvention<TOutput> convention)
        {
            _manager.BeginCheckpoint(key);
            TOutput output;
            string checkpoint = GetCheckpointName();
            string checkpointNumber = GetCheckpointNumber();
            bool hasBeenRunBefore = HasExistingSerializedOutput();
            CheckpointStatus status = _manager.GetCheckpointStatus(hasBeenRunBefore);
            if (status.ShouldRun())
            {
                IDirectoryLocation tempDir = GetTempDirectory(key, status == CheckpointStatus.Run);
                _logger.Info("Running checkpoint {0}: {1}", checkpointNumber, checkpoint);
                Benchmark benchmark = new Benchmark();
                output = run(input, tempDir);
                _logger.Info("Elapsed time (step/time(sec)/name)\t{0}\t{1:F1}\t{2}",
                    checkpointNumber, benchmark.GetElapsedTime(), checkpoint);
                output = convention.Move(output, (source, destination) => source.MoveAndLink(destination));
                Save(output);
                if (!_retainTemps) tempDir.Delete();
            }
            else
            {
                output = Load<TOutput>();
            }
            _manager.EndCheckpoint();
            return output;
        }

        public TOutput RunCheckpoint<TOutput>(string checkpointName, Func<ICheckpointRunner, IDirectoryLocation, TOutput> run)
        {
            return RunCheckpoint(checkpointName, dir => run(this, dir));
        }

        public TOutput RunCheckpoint<TInput, TOutput>(
            string checkpointName,
            Func<ICheckpointRunner, TInput, IDirectoryLocation, TOutput> run,
            TInput input,
            ILoadingConvention<TInput, TOutput> convention)
        {
            return RunCheckpoint(checkpointName, (runInput, dir) => run(this, runInput, dir), input, convention);
        }

        public TOutput RunCheckpoint<TInput, TOutput>(
            string key,
            Func<TInput, IDirectoryLocation, TOutput> run,
            TInput input,
            ILoadingConvention<TInput, TOutput> convention)
        {
            _manager.BeginCheckpoint(key);
            TOutput output;
            string checkpoint = GetCheckpointName();
            string checkpointNumber = GetCheckpointNumber();
            bool hasBeenRunBefore = HasExistingSerializedOutput();
            CheckpointStatus status = _manager.GetCheckpointStatus(hasBeenRunBefore);
            if (status.ShouldRun())
            {
                IDirectoryLocation tempDir = GetTempDirectory(key, status == CheckpointStatus.Run);
                _logger.Info("Running checkpoint {0}: {1}", checkpointNumber, checkpoint);
                Benchmark benchmark = new Benchmark();
                output = run(input, tempDir);
                _logger.Info("Elapsed time (step/time(sec)/name)\t{0}\t{1:F1}\t{2}",
                    checkpointNumber, benchmark.GetElapsedTime(), checkpoint);
                convention.Move(output, (source, destination) => source.MoveAndLink(destination));
                output = convention.Load(input);
                Save(output);
                if (!_retainTemps) tempDir.Delete();
            }
            else if (!hasBeenRunBefore)
            {
                _logger.Info("Skipping checkpoint {0}: {1} by loading expected results", checkpointNumber, checkpoint);
                output = convention.Load(input);
                Save(output);
            }
            else
            {
                output = Load<TOutput>();
            }
            _manager.EndCheckpoint();
            return output;
        }

        public TOut RunCheckpoint<TOut>(string key, Func<IDirectoryLocation, TOut> function)
        {
            _manager.BeginCheckpoint(key);
            TOut output;
            bool hasExistingSerializedOutput = HasExistingSerializedOutput();
            CheckpointStatus status = _manager.GetCheckpointStatus(hasExistingSerializedOutput);
            if (status.ShouldRun())
            {
                string checkpoint = GetCheckpointName();
                string checkpointNumber = GetCheckpointNumber();
                IDirectoryLocation tempDir = GetTempDirectory(checkpoint, status == CheckpointStatus.Run);
                _logger.Info("Running checkpoint {0}: {1}", checkpointNumber, checkpoint);
                Benchmark benchmark = new Benchmark();
                output = function(tempDir);
                _logger.Info("Elapsed time (step/time(sec)/name)\t{0}\t{1:F1}\t{2}",
                    checkpointNumber, benchmark.GetElapsedTime(), checkpoint);
                Save(output);
                if (!_retainTemps) tempDir.Delete();
            }
            else
            {
                output = Load<TOut>();
            }
            _manager.EndCheckpoint();
            return output;
        }

        private void Save<T>(T output)
        {
            _serializer.Save(GetSerializedCheckpointName(), output);
        }

        private T Load<T>()
        {
            return _serializer.Load<T>(GetSerializedCheckpointName());
        }

        private bool HasExistingSerializedOutput()
        {
            return _serializer.Exists(GetSerializedCheckpointName());
        }

        private string GetCheckpointName()
        {
            return string.Join(".", _manager.CheckpointNames);
        }

        private string GetCheckpointNumber()
        {
            return string.Join(".", _manager.CheckpointNumbers);
        }

        private string GetZeroPaddedCheckpointNumber()
        {
            return string.Join(".", _manager.CheckpointNumbers.Select(num => string.Format("{0:00}", num)));
        }

        private string GetSerializedCheckpointName()
        {
            return string.Format("{0}-{1}", GetZeroPaddedCheckpointNumber(), GetCheckpointName().RemoveWhiteSpace());
        }

        private IDirectoryLocation GetTempDirectory(string key, bool makeClean)
        {
            var tempDirectory = _placeForTempFolders.GetDirectoryLocation(key.RemoveWhiteSpace() + "Temp");
            if (makeClean)
                tempDirectory.CreateClean();
            return tempDirectory;
        }

        public void Dispose()
        {
            if (_disposed) return;
            _manager.Dispose();
            _disposed = true;
        }
    }
}
