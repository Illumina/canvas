using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using Isas.Shared.Utilities;

namespace Illumina.SecondaryAnalysis.Workflow
{
    /// <summary>
    /// throws StopCheckpointFoundException when the specified stop checkpoint has finished running
    /// </summary>
    public interface ICheckpointManager : IDisposable
    {
        void BeginCheckpoint(string key);
        bool ShouldRun(bool hasRunBefore);
        List<string> CheckpointNames { get; }
        List<int> CheckpointNumbers { get; }
        void EndCheckpoint();
    }

    [Serializable]
    public class StopCheckpointFoundException : Exception
    {
        public StopCheckpointFoundException(string message) : base(message)
        {
        }
        protected StopCheckpointFoundException(SerializationInfo info, StreamingContext ctxt) : base(info, ctxt)
        {
        }
    }

    [Serializable]
    public class StopCheckpointNotFoundException : Exception
    {
        public StopCheckpointNotFoundException(string message) : base(message)
        {
        }
        protected StopCheckpointNotFoundException(SerializationInfo info, StreamingContext ctxt) : base(info, ctxt)
        {
        }
    }

    [Serializable]
    public class StartCheckpointNotFoundException : Exception
    {
        public StartCheckpointNotFoundException(string message) : base(message)
        {
        }
        protected StartCheckpointNotFoundException(SerializationInfo info, StreamingContext ctxt) : base(info, ctxt)
        {
        }
    }

    public class CheckpointManager : ICheckpointManager
    {
        private readonly ILogger _logger;
        //each entry can be either a checkpoint name or checkpoint number
        private readonly List<string> _startCheckpoint;
        private bool _startCheckpointFound;

        //each entry can be either a checkpoint name or checkpoint number
        private readonly List<string> _stopCheckpoint;
        private bool _stopCheckpointFound;

        private readonly List<string> _checkpointNames;
        public List<string> CheckpointNames
        {
            get
            {
                if (_disposed) throw new ObjectDisposedException(GetType().FullName);
                return _checkpointNames;
            }
        }

        private readonly List<int> _checkpointNumbers;
        public List<int> CheckpointNumbers
        {
            get
            {
                if (_disposed) throw new ObjectDisposedException(GetType().FullName);
                return _checkpointNumbers;
            }
        }

        private int _previousCheckpointNumber;
        private string _lastCheckpointComparedToStartingCheckpoint;
        private bool _disposed;

        public CheckpointManager(ILogger logger, string startCheckpoint = null, string stopCheckpoint = null)
        {
            _logger = logger;
            _startCheckpoint = string.IsNullOrEmpty(startCheckpoint) ? new List<string>() : startCheckpoint.Split('.').ToList();
            if (_startCheckpoint.Any(string.IsNullOrWhiteSpace))
                throw new ArgumentException("Starting checkpoints cannot be empty");
            var invalidNumberedCheckpoints = _startCheckpoint.Where(checkpointString =>
            {
                int checkpointNumber;
                if (int.TryParse(checkpointString, out checkpointNumber) && checkpointNumber <= 0)
                    return true;
                return false;
            });
            if (invalidNumberedCheckpoints.Any())
                throw new ArgumentException("Numbered starting checkpoints must be greater than or equal to 1");


            _stopCheckpoint = string.IsNullOrEmpty(stopCheckpoint) ? new List<string>() : stopCheckpoint.Split('.').ToList();
            if (_stopCheckpoint.Any(string.IsNullOrWhiteSpace))
                throw new ArgumentException("Stopping checkpoints cannot be empty");
            invalidNumberedCheckpoints = _stopCheckpoint.Where(checkpointString =>
            {
                int checkpointNumber;
                if (int.TryParse(checkpointString, out checkpointNumber) && checkpointNumber <= 0)
                    return true;
                return false;
            });
            if (invalidNumberedCheckpoints.Any())
                throw new ArgumentException("Numbered stopping checkpoints must be greater than or equal to 1");

            //todo: validate that any numbered stop checkpoint is at or after any numbered start checkpoints given that they have the same parent checkpoint

            _startCheckpointFound = _startCheckpoint.IsNullOrEmpty();
            _stopCheckpointFound = _stopCheckpoint.IsNullOrEmpty();
            _previousCheckpointNumber = 0;
            _checkpointNames = new List<string>();
            _checkpointNumbers = new List<int>();
        }

        public void BeginCheckpoint(string key)
        {
            if (_disposed) throw new ObjectDisposedException(GetType().FullName);
            if (key.IsNullOrWhiteSpace())
            {
                _disposed = true;
                throw new ArgumentException("Checkpoint name cannot be null or whitespace");
            }
            if (key.IsInt())
            {
                _disposed = true;
                throw new ArgumentException(string.Format("Invalid checkpoint name '{0}'. Checkpoint name cannot be an integer", key));
            }
            CheckpointNames.Add(key);
            CheckpointNumbers.Add(_previousCheckpointNumber + 1);
            // we are potentially starting a new level of nested checkpoints. We want the new level to start from checkpoint #1
            _previousCheckpointNumber = 0;

            if (!_startCheckpointFound)
            {
                int checkpointLevel = CheckpointNames.Count;
                _lastCheckpointComparedToStartingCheckpoint = GetCheckpointComparedToStartingCheckpoint();
                if (CurrentCheckpointComparedToStartingCheckpoint.Last().Equals(_startCheckpoint[checkpointLevel - 1], StringComparison.OrdinalIgnoreCase)
                    && checkpointLevel == _startCheckpoint.Count)
                {
                    _logger.Info("Found specified start checkpoint '{0}'", GetStartingCheckpoint());
                    _startCheckpointFound = true;
                }
            }

            if (IsStopCheckpoint && !_startCheckpointFound)
            {
                _disposed = true;
                throw new StartCheckpointNotFoundException(string.Format(
                    "Found stop checkpoint '{0}', but the specified start checkpoint '{1}' was never reached. The specified stop checkpoint must occur after the specified start checkpoint.",
                    GetStopCheckpoint(), GetStartingCheckpoint()));
            }
        }

        public void EndCheckpoint()
        {
            if (_disposed) throw new ObjectDisposedException(GetType().FullName);
            if (IsStopCheckpoint)
            {
                _disposed = true;
                _stopCheckpointFound = true;
                string message = string.Format("Stop checkpoint '{0}' finished.", GetStopCheckpoint());
                _logger.Info(message);
                // the caller must handle this as appropriate
                throw new StopCheckpointFoundException(message);
            }
            int checkpointLevel = CheckpointNames.Count();
            if (!_startCheckpointFound && checkpointLevel != _startCheckpoint.Count())
            {
                _disposed = true;
                throw new StartCheckpointNotFoundException(
                    string.Format("Execution of checkpoint '{0}' finished before finding the specified nested starting checkpoint '{1}'. The last nested checkpoint found was \"{2}\"",
                        GetCheckpointComparedToStartingCheckpoint(), GetStartingCheckpoint(), _lastCheckpointComparedToStartingCheckpoint));
            }
            CheckpointNames.RemoveLast();
            //we potentially finished a level of nested checkpoints, so we reload our checkpoint counter 
            _previousCheckpointNumber = CheckpointNumbers.RemoveLast();
        }

        public bool ShouldRun(bool hasRunBefore)
        {
            // We want to start running code (as opposed to loading from the repository) when we find the 
            //  checkpoint provided during initialization.
            // When the key matches the starting checkpoint, we change state to indicate that we have 
            //  found the starting checkpoint
            if (_startCheckpointFound)
            {
                // if there is no start checkpoint and this step has already been run before we don't want to run it again
                if (_startCheckpoint.IsNullOrEmpty() && hasRunBefore)
                    return false;

                // Once we've found the starting checkpoint, we run all subsequent checkpoints up to the stopping checkpoint
                return true;
            }
            // We haven't found the checkpoint yet, so check if the current checkpoint is part of the starting checkpoint
            int checkpointLevel = CheckpointNames.Count;
            _lastCheckpointComparedToStartingCheckpoint = GetCheckpointComparedToStartingCheckpoint();
            if (CurrentCheckpointComparedToStartingCheckpoint.Last() == _startCheckpoint[checkpointLevel - 1])
            {
                // We've found part of the current checkpoint!
                return true;
            }
            // We still haven't found the starting checkpoint, so don't run the code (load instead).
            return false;
        }

        /// <summary>
        /// each entry will be from one of the corresponding entries in either _currentCheckpointNames or _currentCheckpointNumbers
        /// depending on the type (name or number) of the corresponding starting checkpoint entry
        /// e.g. given:
        ///        _currentCheckpointNames = { A, A2, A2sub3} 
        ///        _currentCheckpointNumbers = {1, 2, 3}
        ///        _startingCheckpoint = {A, 2, A2sub4}
        ///      then:
        ///        _currentCheckpointComparedToStartingCheckpoint = {A2sub3, 2, A}
        /// </summary>
        /// <returns>List of current checkpoint level entries</returns>
        private List<string> CurrentCheckpointComparedToStartingCheckpoint
        {
            get
            {
                List<string> currentCheckpointComparedToStartingCheckpoint = new List<string>();
                for (int checkpointLevelIndex = 0; checkpointLevelIndex < CheckpointNames.Count; checkpointLevelIndex++)
                {
                    int checkpointNumber = CheckpointNumbers[checkpointLevelIndex];
                    string checkpointName = CheckpointNames[checkpointLevelIndex];
                    if (_startCheckpoint[checkpointLevelIndex].IsInt())
                    {
                        currentCheckpointComparedToStartingCheckpoint.Add(checkpointNumber.ToString());
                    }
                    else
                    {
                        currentCheckpointComparedToStartingCheckpoint.Add(checkpointName);
                    }
                }
                return currentCheckpointComparedToStartingCheckpoint;
            }
        }

        private bool IsStopCheckpoint
        {
            get
            {
                for (int checkpointLevelIndex = 0; checkpointLevelIndex < CheckpointNames.Count; checkpointLevelIndex++)
                {
                    if (checkpointLevelIndex >= _stopCheckpoint.Count) return false;
                    int checkpointNumber = CheckpointNumbers[checkpointLevelIndex];
                    string checkpointName = CheckpointNames[checkpointLevelIndex];
                    if (!_stopCheckpoint[checkpointLevelIndex].Equals(checkpointName, StringComparison.OrdinalIgnoreCase) &&
                        _stopCheckpoint[checkpointLevelIndex] != checkpointNumber.ToString())
                    {
                        return false;
                    }
                }
                return true;
            }
        }

        private string GetCheckpointComparedToStartingCheckpoint()
        {
            return string.Join(".", CurrentCheckpointComparedToStartingCheckpoint);
        }

        private string GetStartingCheckpoint()
        {
            return string.Join(".", _startCheckpoint);
        }

        private string GetStopCheckpoint()
        {
            return string.Join(".", _stopCheckpoint);
        }

        /// <summary>
        /// This should be called after all checkpoints are complete
        /// We validate that any specified start/stop checkpoints were actually found
        /// </summary>
        public void Dispose()
        {
            if (_disposed) return;
            if (_startCheckpoint.Any() && !_startCheckpointFound)
                throw new StartCheckpointNotFoundException($"All checkpoints are finished and we failed to find the specified start checkpoint '{GetStartingCheckpoint()}'");
            if (_stopCheckpoint.Any() && !_stopCheckpointFound)
                throw new StopCheckpointNotFoundException($"All checkpoints are finished and we failed to find the specified stop checkpoint '{GetStopCheckpoint()}'");
            _disposed = true;
        }
    }
}

