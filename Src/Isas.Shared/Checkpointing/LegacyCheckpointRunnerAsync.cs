using System;
using System.Diagnostics;
using System.Threading.Tasks;

namespace Isas.Shared.Checkpointing
{
    public static class LegacyCheckpointRunnerAsync
    {
        /// <summary>
        /// Enable legacy checkpoints to use the checkpointer
        /// Legacy checkpoints do not return any results since they assume downstream checkpoints are aware of the hard-coded paths on the filesystem.
        /// Here we assign them a result which is just a string indicating the checkpoint has completed without returning any output
        /// </summary>
        /// <param name="runner"></param>
        /// <param name="name"></param>
        /// <param name="a"></param>
        public static Task RunCheckpointAsync(this ICheckpointRunnerAsync runner, string name, Action a)
        {
            string result = $"{name} complete. No output from this checkpoint";
            return runner.RunCheckpointAsync(name, a, result);
        }

        public static Task<T> RunCheckpointAsync<T>(this ICheckpointRunnerAsync runner, string name, Action a, T result) where T : class
        {
            Func<T> func = () =>
            {
                a.Invoke();
                return result;
            };
            ISandboxedCheckpoint<string, T> sandboxedCheckpoint = new LegacyCheckpointRunner.NullSandboxedCheckpoint<string, T>(func);
            return runner.RunCheckpointAsync(name, sandboxedCheckpoint, null, new LegacyCheckpointRunner.NullLoadingConvention<string, T>(result));
        }
    }
}