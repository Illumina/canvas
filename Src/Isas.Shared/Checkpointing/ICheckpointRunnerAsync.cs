using System;
using System.Threading.Tasks;

namespace Isas.Shared.Checkpointing
{
    public interface ICheckpointRunnerAsync : ICheckpointRunner
    {
        Task<TResult> RunCheckpointAsync<TResult>(string checkpointName, Func<ICheckpointRunnerAsync, IDirectoryLocation, TResult> run);
        Task<TResult> RunCheckpointAsync<TResult>(string checkpointName, Func<TResult> func);
        Task<TResult> RunCheckpointAsync<T, TResult>(string checkpointName, Func<T, IDirectoryLocation, TResult> run, T input, INamingConvention<TResult> convention);
        Task<TResult> RunCheckpointAsync<T, TResult>(string checkpointName, Func<T, IDirectoryLocation, TResult> run, T input, ILoadingConvention<T, TResult> convention);
    }

    public static class CheckpointRunnerAsyncExtensions
    {
        public static Task<TOutput> RunCheckpointAsync<TInput, TOutput>(this ICheckpointRunnerAsync runner,
    string key,
    ISandboxedCheckpoint<TInput, TOutput> wrapper,
    TInput input,
    ILoadingConvention<TInput, TOutput> loadingConvention)
        {
            return runner.RunCheckpointAsync(key, wrapper.Run, input, loadingConvention);
        }
        public static Task<TOutput> RunCheckpointAsync<TConvention, TInput, TOutput>(this ICheckpointRunnerAsync runner,
    ILoadableResultCheckpoint<TConvention, TInput, TOutput> wrapper,
    TInput input,
    TConvention convention)
        {
            return runner.RunCheckpointAsync(wrapper, input, wrapper.GetLoadingConvention(convention));
        }

        public static Task<TOutput> RunCheckpointAsync<TInput, TOutput>(this ICheckpointRunnerAsync runner,
    IAbstractCheckpoint<TInput, TOutput> wrapper,
    TInput input,
    ILoadingConvention<TInput, TOutput> loadingConvention)
        {
            return runner.RunCheckpointAsync(wrapper.GetStepName(), wrapper, input, loadingConvention);
        }
    }
}