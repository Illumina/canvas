using System;
using System.Threading.Tasks;

namespace Isas.Shared.Checkpointing
{
    public interface ICheckpointRunner : IDisposable
    {
        /// <summary>
        /// Run a checkpoint
        ///  Either the <paramref name="function"/> is run, and the result 
        ///  of type <typeparamref name="TResult"/> is saved and returned;
        ///  or the result is loaded.
        /// </summary>
        /// <typeparam name="TResult"></typeparam>
        /// <param name="key">The name of the checkpoint being run.</param>
        /// <param name="function">The code comprising the checkpoint. 
        ///  Must return an object of type <typeparamref name="TResult"/>. 
        ///  This function should execute all desired steps (such as moving
        ///  of results to final locations) before the checkpoint is saved.</param>
        /// <returns>The result of executing <paramref name="function"/></returns>
        TResult RunCheckpoint<TResult>(string key, Func<TResult> function);

        /// <summary>
        /// Run a checkpoint in a temporary directory.
        /// Either the <paramref name="function"/> is run, and the result 
        ///  of type <typeparamref name="TResult"/> is saved and returned;
        ///  or the result is loaded.
        /// </summary>
        /// <typeparam name="TResult"></typeparam>
        /// <param name="key">The name of the checkpoint being run.</param>
        /// <param name="function">The code comprising the checkpoint. 
        ///  Must return an object of type <typeparamref name="TResult"/>. 
        ///  This function should execute all desired steps (such as moving
        ///  of results to final locations) before the checkpoint is saved.
        /// This function accepts an IDirectoryLocation as input - this directory
        /// is a temporary directory managed by the checkpoint manager.</param>
        /// <returns>The result of executing <paramref name="function"/></returns>
        TResult RunCheckpoint<TResult>(string key, Func<IDirectoryLocation, TResult> function);
        TResult RunCheckpoint<TResult>(string checkpointName, Func<ICheckpointRunner, IDirectoryLocation, TResult> run);
        TResult RunCheckpoint<T, TResult>(string key,
            Func<T, IDirectoryLocation, TResult> run,
            T input,
            ILoadingConvention<T, TResult> loadingConvention);
        TResult RunCheckpoint<T, TResult>(string key,
            Func<T, IDirectoryLocation, TResult> run,
            T input,
            INamingConvention<TResult> namingConvention);

        TOut RunCheckpoint<TIn, TOut>(string key, Func<ICheckpointRunner, TIn, IDirectoryLocation, TOut> run,
            TIn input, ILoadingConvention<TIn, TOut> loadingConvention);
    }

    public interface ICheckpointRunnerAsync : ICheckpointRunner
    {
        Task<TResult> RunCheckpointAsync<TResult>(string checkpointName, Func<ICheckpointRunnerAsync, IDirectoryLocation, TResult> run);
        Task<TResult> RunCheckpointAsync<TResult>(string checkpointName, Func<TResult> func);
        Task<TResult> RunCheckpointAsync<T, TResult>(string checkpointName, Func<T, IDirectoryLocation, TResult> run, T input, INamingConvention<TResult> convention);
        Task<TResult> RunCheckpointAsync<T, TResult>(string checkpointName, Func<T, IDirectoryLocation, TResult> run, T input, ILoadingConvention<T, TResult> convention);
    }

    /// <summary>
    /// An interface for loading results from disk. This interface is used by the ISkippableCheckpoint interface
    /// The following is a requirement for any class implementing this interface:
    ///       var result = Move(source, (source, destination) => source.MoveTo(destination));
    ///       var loadedResult = Load();
    ///       Assert.Equal(result, loadedResult)
    /// This can be achieved by having the Move implementation directly return the result of invoking the Load method
    /// Ultimately however, whatever files are moved in the Move method should also be loaded by the Load method
    /// </summary>
    /// <typeparam name="TResult"></typeparam>
    /// <typeparam name="TInput"></typeparam>
    public interface ILoadingConvention<in TInput, TResult>
    {
        TResult Load(TInput input);
        void Move(TResult source, Action<IFileLocation, IFileLocation> move);
    }

    public interface INamingConvention<TResult>
    {
        TResult Move(TResult source, Action<IFileLocation, IFileLocation> move);
    }

    public interface IMoveableResultCheckpoint<in T, in TInput, TResult> : IAbstractCheckpoint<TInput, TResult>
    {
        INamingConvention<TResult> GetNamingConvention(T fileNameConventionSetting);
    }

    public interface ILoadableResultCheckpoint<in T, in TInput, TResult> : IAbstractCheckpoint<TInput, TResult>
    {
        ILoadingConvention<TInput, TResult> GetLoadingConvention(T fileNameConventionSetting);
    }

    public interface IAbstractCheckpoint<in TInput, out TOutput> : ISandboxedCheckpoint<TInput, TOutput>
    {
        string GetStepName();
    }

    public interface ISandboxedCheckpoint<in TInput, out TOutput>
    {
        TOutput Run(TInput input, IDirectoryLocation sandbox);
    }

    // Delegates for extension methods

    public delegate void OneOut<Tout>(out Tout output);

    public delegate void NestedOneOut<Tout>(ICheckpointRunner runner, out Tout output);

    public delegate void TwoOuts<T1, T2>(out T1 output1, out T2 output2);

    public delegate void NestedTwoOuts<T1, T2>(ICheckpointRunner runner, out T1 output1, out T2 output2);

    public static class CheckpointRunnerExtensions
    {



        public static void RunCheckpoint(this ICheckpointRunner runner, string name, Action<ICheckpointRunner> a)
        {
            string result = $"{name} complete. No output from this checkpoint";
            Func<ICheckpointRunner, IDirectoryLocation, string> func = (checkpointer, dir) =>
            {
                a.Invoke(checkpointer);
                return result;
            };
            runner.RunCheckpoint(name, func);
        }

        public static TOutput RunCheckpoint<TInput, TOutput>(this ICheckpointRunner runner,
            string key,
            ISandboxedCheckpoint<TInput, TOutput> wrapper,
            TInput input,
            INamingConvention<TOutput> namingConvention)
        {
            return runner.RunCheckpoint(key, wrapper.Run, input, namingConvention);
        }

        public static TOutput RunCheckpoint<TInput, TOutput>(this ICheckpointRunner runner,
            string key,
            ISandboxedCheckpoint<TInput, TOutput> wrapper,
            TInput input,
            ILoadingConvention<TInput, TOutput> loadingConvention)
        {
            return runner.RunCheckpoint(key, wrapper.Run, input, loadingConvention);
        }

        public static TOutput RunCheckpoint<TInput, TOutput>(this ICheckpointRunner runner,
            IAbstractCheckpoint<TInput, TOutput> wrapper,
            TInput input,
            INamingConvention<TOutput> namingConvention)
        {
            return runner.RunCheckpoint(wrapper.GetStepName(), wrapper, input, namingConvention);
        }

        public static TOutput RunCheckpoint<TInput, TOutput>(this ICheckpointRunner runner,
            IAbstractCheckpoint<TInput, TOutput> wrapper,
            TInput input,
            ILoadingConvention<TInput, TOutput> loadingConvention)
        {
            return runner.RunCheckpoint(wrapper.GetStepName(), wrapper, input, loadingConvention);
        }

        public static TOutput RunCheckpoint<TConvention, TInput, TOutput>(this ICheckpointRunner runner,
            IMoveableResultCheckpoint<TConvention, TInput, TOutput> wrapper,
            TInput input,
            TConvention convention)
        {
            return runner.RunCheckpoint(wrapper, input, wrapper.GetNamingConvention(convention));
        }

        public static TOutput RunCheckpoint<TConvention, TInput, TOutput>(this ICheckpointRunner runner,
            ILoadableResultCheckpoint<TConvention, TInput, TOutput> wrapper,
            TInput input,
            TConvention convention)
        {
            return runner.RunCheckpoint(wrapper, input, wrapper.GetLoadingConvention(convention));
        }

        public static void RunCheckpoint<TOut>(this ICheckpointRunner runner,
            string key, out TOut output, OneOut<TOut> function)
        {
            output = runner.RunCheckpoint(key, () =>
            {
                TOut out1;
                function(out out1);
                return out1;
            });
        }

        public static void RunCheckpoint<TOut>(this ICheckpointRunner runner,
            string key, out TOut output, NestedOneOut<TOut> function)
        {
            output = runner.RunCheckpoint(key, (checkpointer, dir) =>
            {
                TOut out1;
                function(checkpointer, out out1);
                return out1;
            });
        }

        public static void RunCheckpoint<TOut1, TOut2>(this ICheckpointRunner runner,
            string key, out TOut1 output1, out TOut2 output2, TwoOuts<TOut1, TOut2> function)
        {
            Tuple<TOut1, TOut2> tup;
            runner.RunCheckpoint(key, out tup, delegate (out Tuple<TOut1, TOut2> _tup)
            {
                TOut1 out1;
                TOut2 out2;
                function(out out1, out out2);
                _tup = new Tuple<TOut1, TOut2>(out1, out2);
            });
            output1 = tup.Item1;
            output2 = tup.Item2;
        }

        public static void RunCheckpoint<TOut1, TOut2>(this ICheckpointRunner runner,
            string key, out TOut1 output1, out TOut2 output2, NestedTwoOuts<TOut1, TOut2> function)
        {
            Tuple<TOut1, TOut2> tup;
            runner.RunCheckpoint(key, out tup, delegate (ICheckpointRunner checkpointer, out Tuple<TOut1, TOut2> _tup)
            {
                TOut1 out1;
                TOut2 out2;
                function(checkpointer, out out1, out out2);
                _tup = new Tuple<TOut1, TOut2>(out1, out2);
            });
            output1 = tup.Item1;
            output2 = tup.Item2;
        }
    }
}