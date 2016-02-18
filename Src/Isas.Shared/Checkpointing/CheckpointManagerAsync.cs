using System;
using Illumina.SecondaryAnalysis;
using ILMNcommon.Common;

namespace Isas.Shared.Checkpointing
{
    public class CheckpointManagerAsync : IDisposable
    {
        private readonly ILogger _logger;
        private readonly string _startCheckpoint;
        private readonly string _stopCheckpoint;
        private bool _startCheckpointFound;
        private int _currentCheckpointNumber;
        private bool _stopCheckpointFound;
        private bool _disposed;

        public CheckpointManagerAsync(ILogger logger, string startCheckpoint = null, string stopCheckpoint = null)
        {
            _logger = logger;
            _startCheckpoint = startCheckpoint;
            _stopCheckpoint = stopCheckpoint;

            string error;
            if (!IsValidCheckpoint(_startCheckpoint, out error))
                throw new ArgumentException($"Invalid start checkpoint '{_startCheckpoint}': {error}");
            if (!IsValidCheckpoint(_stopCheckpoint, out error))
                throw new ArgumentException($"Invalid stop checkpoint '{_stopCheckpoint}': {error}");
        }

        private static bool IsValidCheckpoint(string checkpoint, out string error)
        {
            error = null;
            if (checkpoint == null) return true;
            if (checkpoint.RemoveWhiteSpace().Empty())
            {
                error = "checkpoint cannot be empty";
                return false;
            }
            int checkpointNumber;
            if (int.TryParse(checkpoint, out checkpointNumber) && checkpointNumber <= 0)
            {
                error = "checkpoint must be greater than or equal to 1";
                return false;
            }
            return true;
        }

        public Checkpoint CreateCheckpoint(string checkpoint)
        {
            if (_disposed) throw new ObjectDisposedException(GetType().FullName);
            if (checkpoint.IsNullOrWhiteSpace())
                throw new ArgumentException("checkpoint name cannot be null or whitespace");

            _currentCheckpointNumber++;
            bool isStartCheckpoint = false;
            if (IsMatchingCheckpoint(checkpoint, _currentCheckpointNumber, _startCheckpoint))
            {
                _logger.Info($"Found specified start checkpoint {_currentCheckpointNumber}: '{checkpoint}'");
                isStartCheckpoint = true;
                _startCheckpointFound = true;
            }
            if (IsMatchingCheckpoint(checkpoint, _currentCheckpointNumber, _stopCheckpoint))
                _stopCheckpointFound = true;

            if (_stopCheckpointFound && _startCheckpoint != null && !_startCheckpointFound)
                throw new StartCheckpointNotFoundException(
                    $"Found stop checkpoint '{_stopCheckpoint}', but the specified start checkpoint '{_startCheckpoint}' was never reached. The specified stop checkpoint must occur after the specified start checkpoint or be the specified start checkpoint itself.");

            bool startCheckpointExists = _startCheckpoint != null;
            bool beforeStartCheckpoint = startCheckpointExists && !_startCheckpointFound;
            return new Checkpoint(_currentCheckpointNumber, checkpoint, beforeStartCheckpoint, startCheckpointExists, isStartCheckpoint, _stopCheckpointFound);
        }

        private static bool IsMatchingCheckpoint(string queryCheckpointName, int queryCheckpointNumber, string checkpointNameOrNumber)
        {
            return queryCheckpointNumber.ToString() == checkpointNameOrNumber || queryCheckpointName.Equals(checkpointNameOrNumber, StringComparison.OrdinalIgnoreCase);
        }

        /// <summary>
        /// This should be called after all checkpoints are complete
        /// We validate that any specified start/stop checkpoints were actually found
        /// </summary>
        public void Dispose()
        {
            if (_disposed) return;
            if (_startCheckpoint != null && !_startCheckpointFound)
                throw new StartCheckpointNotFoundException($"All checkpoints are finished and we failed to find the specified start checkpoint '{_startCheckpoint}'");
            if (_stopCheckpoint != null && !_stopCheckpointFound)
                throw new StopCheckpointNotFoundException($"All checkpoints are finished and we failed to find the specified stop checkpoint '{_stopCheckpoint}'");
            _disposed = true;
        }

        public CheckpointManagerAsync CreateChildCheckpointerManager(string startCheckpoint, string stopCheckpoint)
        {
            if (_disposed) throw new ObjectDisposedException(GetType().FullName);
            return new CheckpointManagerAsync(_logger, startCheckpoint, stopCheckpoint);
        }
    }
}