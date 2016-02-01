using System;

namespace Isas.Shared.Checkpointing
{
    /// <summary>
    /// This class contains all the information about a single checkpoint such as:
    ///     -checkpoint number
    ///     -checkpoint name
    ///     -whether it is the start checkpoint
    ///     -whether it is the stop checkpoint
    ///     -whether the step should run or results should be loaded from the json checkpoint file.
    /// The "Start checkpoint" is passed from the Isas command-line using the -c option
    /// The "Stop checkpoint" is passed from the Isas command-line using the -s option
    /// </summary>
    public class Checkpoint
    {
        private readonly bool _beforeStartCheckpoint;
        private readonly bool _startCheckpointExists;
        public int Number { get; }
        public string Name { get; }
        public bool IsStartCheckpoint { get; set; }

        public bool ShouldRun(bool hasBeenRunBefore)
        {
            // haven't found the start checkpoint yet
            if (_beforeStartCheckpoint) return false;

            // no start checkpoint specified, so we resume after the last checkpoint that exists
            if (!_startCheckpointExists && hasBeenRunBefore)
                return false;

            return true;
        }

        public bool IsStopCheckpoint { get; }

        public Checkpoint(int number, string name, bool beforeStartCheckpoint, bool startCheckpointExists, bool isStartCheckpoint, bool isStopCheckpoint)
        {
            if (beforeStartCheckpoint && !startCheckpointExists)
                throw new ArgumentException("checkpoint cannot be before the start checkpoint if there is no start checkpoint");
            if (isStartCheckpoint && !startCheckpointExists)
                throw new ArgumentException("checkpoint cannot be the start checkpoint if there is no start checkpoint");
            if (beforeStartCheckpoint && isStartCheckpoint)
                throw new ArgumentException("checkpoint cannot be the start checkpoint if it is before the start checkpoint");
            _beforeStartCheckpoint = beforeStartCheckpoint;
            _startCheckpointExists = startCheckpointExists;
            Number = number;
            Name = name;
            IsStartCheckpoint = isStartCheckpoint;
            IsStopCheckpoint = isStopCheckpoint;
        }

        public override string ToString()
        {
            return $"{Number}: '{Name}'";
        }
    }
}