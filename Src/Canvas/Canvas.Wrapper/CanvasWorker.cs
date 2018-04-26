using System;
using System.Threading.Tasks;
using Illumina.Common.FileSystem;
using Isas.Framework.Checkpointing;
using Isas.Framework.Checkpointing.Legacy;
using Isas.Framework.DataTypes;

namespace Canvas.Wrapper
{
    public delegate IFileLocation SampleStubNamingConvention(SampleInfo sample);

    [System.Obsolete("ICanvasWorker needs to be refactored")]
    public interface ICanvasWorker<TCanvasInput, TCanvasOutput> where TCanvasOutput : ICanvasOutput
    {
        SampleSet<TCanvasOutput> Run(SampleSet<TCanvasInput> inputs, SampleStubNamingConvention sampleStubNamingConvention, ICheckpointRunner checkpointRunner);
        Task<SampleSet<TCanvasOutput>> RunAsync(SampleSet<TCanvasInput> inputs, SampleStubNamingConvention sampleStubNamingConvention, ICheckpointRunner checkpointRunner);
    }

    [Obsolete("ICanvasWorker needs to be refactored")]
    public class CanvasWorker<TCanvasInput, TCanvasOutput> : ICanvasWorker<TCanvasInput, TCanvasOutput> where TCanvasInput : ICanvasCheckpointInput where TCanvasOutput : ICanvasOutput
    {
        private readonly ICanvasCheckpoint<TCanvasInput, TCanvasOutput> _canvasCheckpoint;
        public CanvasWorker(ICanvasCheckpoint<TCanvasInput, TCanvasOutput> canvasCheckpoint)
        {
            _canvasCheckpoint = canvasCheckpoint;
        }

        public SampleSet<TCanvasOutput> Run(SampleSet<TCanvasInput> inputs, SampleStubNamingConvention sampleStubNamingConvention, ICheckpointRunner checkpointRunner)
        {
            return checkpointRunner.RunCheckpoint(_canvasCheckpoint, inputs, sampleStubNamingConvention);
        }

        public Task<SampleSet<TCanvasOutput>> RunAsync(SampleSet<TCanvasInput> inputs, SampleStubNamingConvention sampleStubNamingConvention, ICheckpointRunner checkpointRunner)
        {
            return checkpointRunner.RunCheckpointAsync(_canvasCheckpoint, inputs, sampleStubNamingConvention);
        }
    }

    [Obsolete("ICanvasWorker needs to be refactored")]
    public class NullCanvasWorker<TCanvasInput, TCanvasOutput> : ICanvasWorker<TCanvasInput, TCanvasOutput> where TCanvasOutput : ICanvasOutput
    {
        public Task<SampleSet<TCanvasOutput>> RunAsync(SampleSet<TCanvasInput> inputs, SampleStubNamingConvention sampleStubNamingConvention, ICheckpointRunner checkpointRunner)
        {
            return Task.FromResult(inputs.SelectData(input => default(TCanvasOutput)));
        }

        public SampleSet<TCanvasOutput> Run(SampleSet<TCanvasInput> inputs, SampleStubNamingConvention sampleStubNamingConvention, ICheckpointRunner checkpointRunner)
        {
            return inputs.SelectData(input => default(TCanvasOutput));
        }
    }

    public interface IAbstractCheckpoint<in TInput, out TOutput> : ISandboxedCheckpoint<TInput, TOutput>
    {
        /// <summary>
        /// Get the step name for this checkpoint
        /// </summary>
        /// <returns></returns>
        string GetStepName();
    }

    public interface ISandboxedCheckpoint<in TInput, out TOutput>
    {
        /// <summary>
        /// Do work given input in a temporary directory to produce output
        /// </summary>
        /// <param name="input"></param>
        /// <param name="sandbox"></param>
        /// <returns></returns>
        TOutput Run(TInput input, IDirectoryLocation sandbox);
    }

    /// <summary>
    /// Support running legacy asynchronous checkpoints using ICheckpointRunner interface
    /// </summary>
    public static class CheckpointRunnerAsyncExtensions
    {
        public static TResult RunCheckpoint<TResult>(this ICheckpointRunner checkpointRunner, string checkpointName, Func<ICheckpointRunner, IDirectoryLocation, IFileMover, TResult> run)
        {
            return checkpointRunner.RunCheckpoint(checkpointName, (tempDir, fileMover) => run(checkpointRunner, tempDir, fileMover));
        }

        public static Task<TResult> RunCheckpointAsync<TResult>(this ICheckpointRunner checkpointRunner, string checkpointName, Func<TResult> function)
        {
            return checkpointRunner.RunCheckpoint(checkpointName, () => Task.Run(function));
        }

        public static Task<TResult> RunCheckpointAsync<TResult>(this ICheckpointRunner checkpointRunner, string key, Func<TResult> run, Func<TResult> load)
        {
            return checkpointRunner.RunCheckpoint(key, () => Task.Run(run), () => Task.Run(load));
        }

        public static Task<TResult> RunCheckpointAsync<TResult>(this ICheckpointRunner checkpointRunner, string checkpointName, Func<IDirectoryLocation, IFileMover, TResult> function)
        {
            return checkpointRunner.RunCheckpoint(checkpointName, (tempDir, fileMover) => Task.Run(() => function(tempDir, fileMover)));
        }

        public static Task<TResult> RunCheckpointAsync<TResult>(this ICheckpointRunner checkpointRunner, string checkpointName, Func<IDirectoryLocation, IFileMover, TResult> function, Func<TResult> load)
        {
            return checkpointRunner.RunCheckpoint(checkpointName, (tempDir, fileMover) => Task.Run(() => function(tempDir, fileMover)), () => Task.Run(load));
        }
    }

    public static class LoadingAndMovingCheckpointerAsync
    {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="runner"></param>
        /// <param name="key"></param>
        /// <param name="wrapper"></param>
        /// <param name="input"></param>
        /// <param name="loadingConvention"></param>
        /// <typeparam name="TInput"></typeparam>
        /// <typeparam name="TOutput"></typeparam>
        /// <returns></returns>
        /// <exception cref="NotImplementedException"></exception>
        public static Task<TOutput> RunCheckpointAsync<TInput, TOutput>(this ICheckpointRunner runner, string key, ISandboxedCheckpoint<TInput, TOutput> wrapper, TInput input, ILoadingConvention<TInput, TOutput> loadingConvention)
        {
            var asyncLoadingConvention = new LoadingConventionAsync<TInput, TOutput>(loadingConvention);
            return runner.RunCheckpoint(key, (tempInput, dir) => Task.Run(() => wrapper.Run(tempInput, dir)), input, asyncLoadingConvention);
        }

        public static Task<TOutput> RunCheckpointAsync<TInput, TOutput>(this ICheckpointRunner runner, string key, ISandboxedCheckpoint<TInput, TOutput> wrapper, TInput input,
    INamingConvention<TOutput> namingConvention)
        {
            return runner.RunCheckpointAsync(key, dir => wrapper.Run(input, dir), namingConvention);
        }

        public static Task<TOutput> RunCheckpointAsync<TOutput>(this ICheckpointRunner runner,
string key,
Func<IDirectoryLocation, TOutput> run,
INamingConvention<TOutput> namingConvention)
        {
            var asyncNamingConvention = new NamingConventionAsync<TOutput>(namingConvention);
            return runner.RunCheckpoint(key, dir => Task.Run(() => run(dir)), asyncNamingConvention);
        }

        public static Task<TOutput> RunCheckpointAsync<TConvention, TInput, TOutput>(this ICheckpointRunner runner,
    ILoadableResultCheckpoint<TConvention, TInput, TOutput> wrapper,
    TInput input,
    TConvention convention)
        {
            return runner.RunCheckpointAsync(wrapper, input, wrapper.GetLoadingConvention(convention));
        }

        public static Task<TOutput> RunCheckpointAsync<TConvention, TInput, TOutput>(this ICheckpointRunner runner,
    IMoveableResultCheckpoint<TConvention, TInput, TOutput> wrapper,
    TInput input,
    TConvention convention)
        {
            return runner.RunCheckpointAsync(wrapper, input, wrapper.GetNamingConvention(convention));
        }

        public static Task<TOutput> RunCheckpointAsync<TInput, TOutput>(this ICheckpointRunner runner,
    IAbstractCheckpoint<TInput, TOutput> wrapper,
    TInput input,
    ILoadingConvention<TInput, TOutput> loadingConvention)
        {
            return runner.RunCheckpointAsync(wrapper.GetStepName(), wrapper, input, loadingConvention);
        }

        public static Task<TOutput> RunCheckpointAsync<TInput, TOutput>(this ICheckpointRunner runner,
IAbstractCheckpoint<TInput, TOutput> wrapper,
TInput input,
INamingConvention<TOutput> namingConvention)
        {
            return runner.RunCheckpointAsync(wrapper.GetStepName(), wrapper, input, namingConvention);
        }


        private class LoadingConventionAsync<TInput, TResult> : ILoadingConvention<TInput, Task<TResult>>
        {
            private readonly ILoadingConvention<TInput, TResult> _loadingConvention;

            public LoadingConventionAsync(ILoadingConvention<TInput, TResult> loadingConvention)
            {
                _loadingConvention = loadingConvention;
            }

            public Task<TResult> Load(TInput input)
            {
                return Task.Run(() => _loadingConvention.Load(input));
            }

            public void Move(Task<TResult> source, Action<IFileLocation, IFileLocation> move)
            {
                _loadingConvention.Move(source.Result, move);
            }
        }

        private class NamingConventionAsync<TResult> : INamingConvention<Task<TResult>>
        {
            private readonly INamingConvention<TResult> _namingConvention;

            public NamingConventionAsync(INamingConvention<TResult> namingConvention)
            {
                _namingConvention = namingConvention;
            }

            /// <summary>
            /// Move the source given the provided move action to generate the result
            /// </summary>
            /// <param name="source"></param>
            /// <param name="move"></param>
            /// <returns></returns>
            public Task<TResult> Move(Task<TResult> source, Action<IFileLocation, IFileLocation> move)
            {
                return Task.Run(() => _namingConvention.Move(source.Result, move));
            }
        }
    }
}
