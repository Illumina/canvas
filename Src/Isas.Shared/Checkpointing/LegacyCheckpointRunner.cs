using System;

namespace Isas.Shared.Checkpointing
{
    public static class LegacyCheckpointRunner
    {
        /// <summary>
        /// Enable legacy checkpoints to use the checkpointer
        /// Legacy checkpoints do not return any results since they assume downstream checkpoints are aware of the hard-coded paths on the filesystem.
        /// Here we assign them a result which is just a string indicating the checkpoint has completed without returning any output
        /// </summary>
        /// <param name="runner"></param>
        /// <param name="name"></param>
        /// <param name="a"></param>
        public static void RunCheckpoint(this ICheckpointRunner runner, string name, Action a)
        {
            string result = $"{name} complete. No output from this checkpoint";
            Func<string> func = () =>
            {
                a.Invoke();
                return result;
            };
            NullSandboxedCheckpoint<string, string> sandboxedCheckpoint = new NullSandboxedCheckpoint<string, string>(func);
            runner.RunCheckpoint(name, sandboxedCheckpoint, null, new NullLoadingConvention<string, string>(result));
        }

        /// <summary>
        /// An adapter to allow running a legacy checkpoint using the sandboxed checkpoint functionality.
        /// The temp directory gets created but is not used
        /// </summary>
        public class NullSandboxedCheckpoint<TInput, TOutput> : ISandboxedCheckpoint<TInput, TOutput>
        {
            private readonly Func<TOutput> _output;

            public NullSandboxedCheckpoint(Func<TOutput> output)
            {
                _output = output;
            }

            public TOutput Run(TInput input, IDirectoryLocation tempDirectory)
            {
                return _output.Invoke();
            }
        }

        /// <summary>
        /// An adapter to allow skipping a legacy checkpoint. 
        /// Since the legacy checkpoint does not return any results the Load method simply return some precomputed output that is passed via the constructor
        /// </summary>
        public class NullLoadingConvention<TInput, TOutput> : ILoadingConvention<TInput, TOutput> where TOutput : class
        {
            private readonly TOutput _result;

            public NullLoadingConvention(TOutput result)
            {
                _result = result;
            }

            public void Move(TOutput source, Action<IFileLocation, IFileLocation> a)
            {
            }

            public TOutput Load(TInput input)
            {
                return _result;
            }
        }
    }
}
