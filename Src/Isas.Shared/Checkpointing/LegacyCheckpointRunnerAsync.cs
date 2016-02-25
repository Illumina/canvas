using System;
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
        public static async Task RunCheckpointAsync(this ICheckpointRunnerAsync runner, string name, Action a)
        {
            string result = $"{name} complete. No output from this checkpoint";
            Func<string> func = () =>
            {
                a.Invoke();
                return result;
            };
            ISandboxedCheckpoint<string, string> sandboxedCheckpoint = new LegacyCheckpointRunner.NullSandboxedCheckpoint<string, string>(func);
            await runner.RunCheckpointAsync(name, sandboxedCheckpoint, null, new LegacyCheckpointRunner.NullLoadingConvention<string, string>(result));
        }
    }
}