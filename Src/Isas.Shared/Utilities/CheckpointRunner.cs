using System;
using System.Linq;
using Isas.Shared;
using Isas.Shared.Utilities;

namespace Illumina.SecondaryAnalysis.Workflow
{
    public interface ICheckpointRunner : IDisposable
    {
        /// <summary>
        /// Run a checkpoint
        ///  Either the <paramref name="function"/> is run, and the result 
        ///  of type <typeparamref name="Tout"/> is saved and returned;
        ///  or the result is loaded.
        /// </summary>
        /// <typeparam name="Tout"></typeparam>
        /// <param name="key">The name of the checkpoint being run.</param>
        /// <param name="function">The code comprising the checkpoint. 
        ///  Must return an object of type <typeparamref name="Tout"/>. 
        ///  This function should execute all desired steps (such as moving
        ///  of results to final locations) before the checkpoint is saved.</param>
        /// <returns>The result of executing <paramref name="function"/></returns>
        Tout RunCheckpoint<Tout>(string key, Func<Tout> function);

        /// <summary>
        /// Run a checkpoint in a temporary directory.
        /// Either the <paramref name="function"/> is run, and the result 
        ///  of type <typeparamref name="Tout"/> is saved and returned;
        ///  or the result is loaded.
        /// </summary>
        /// <typeparam name="Tout"></typeparam>
        /// <param name="key">The name of the checkpoint being run.</param>
        /// <param name="function">The code comprising the checkpoint. 
        ///  Must return an object of type <typeparamref name="Tout"/>. 
        ///  This function should execute all desired steps (such as moving
        ///  of results to final locations) before the checkpoint is saved.
        /// This function accepts an IDirectoryLocation as input - this directory
        /// is a temporary directory managed by the checkpoint manager.</param>
        /// <returns>The result of executing <paramref name="function"/></returns>
        Tout RunCheckpoint<Tout>(string key, Func<IDirectoryLocation, Tout> function);

        Tout RunCheckpoint<Tin, Tout>(string key,
            ISandboxedCheckpoint<Tin, Tout> wrapper,
            Tin input,
            ILoadingConvention<Tin, Tout> loadingConvention);
        Tout RunCheckpoint<Tin, Tout>(string key,
            ISandboxedCheckpoint<Tin, Tout> wrapper,
            Tin input,
            INamingConvention<Tout> namingConvention);
    }

    public static class CheckpointRunnerExtensions
    {
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
    }

    public static class ICheckpointRunnerExtenstions
    {
        public static void RunCheckpoint<Tout>(this ICheckpointRunner runner,
            string key, out Tout output, OneOut<Tout> function)
        {
            output = runner.RunCheckpoint<Tout>(key, () =>
            {
                Tout out1;
                function(out out1);
                return out1;
            });
        }

        public static void RunCheckpoint<Tout1, Tout2>(this ICheckpointRunner runner,
            string key, out Tout1 output1, out Tout2 output2, TwoOuts<Tout1, Tout2> function)
        {
            Tuple<Tout1, Tout2> tup;
            runner.RunCheckpoint(key, out tup, delegate (out Tuple<Tout1, Tout2> _tup)
            {
                Tout1 out1;
                Tout2 out2;
                function(out out1, out out2);
                _tup = new Tuple<Tout1, Tout2>(out1, out2);
            });
            output1 = tup.Item1;
            output2 = tup.Item2;
        }
    }

    public interface ISandboxedCheckpoint<TInput, TOutput>
    {
        TOutput Run(TInput input, IDirectoryLocation sandbox);
    }

    public interface IAbstractCheckpoint<TInput, TOutput> : ISandboxedCheckpoint<TInput, TOutput>
    {
        string GetStepName();
    }

    public interface ILoadableResultCheckpoint<T, TInput, TResult> : IAbstractCheckpoint<TInput, TResult>
    {
        ILoadingConvention<TInput, TResult> GetLoadingConvention(T fileNameConventionSetting);
    }

    public interface IMoveableResultCheckpoint<T, TInput, TResult> : IAbstractCheckpoint<TInput, TResult>
    {
        INamingConvention<TResult> GetNamingConvention(T fileNameConventionSetting);
    }

    public interface INamingConvention<TResult>
    {
        TResult Move(TResult source, Action<IFileLocation, IFileLocation> move);
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
    public interface ILoadingConvention<TInput, TResult>
    {
        TResult Load(TInput input);
        void Move(TResult source, Action<IFileLocation, IFileLocation> move);
    }

    // Delegates for extension methods
    public delegate void OneOut<Tout>(out Tout output);

    public delegate void TwoOuts<T1, T2>(out T1 output1, out T2 output2);

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

        public Tout RunCheckpoint<Tout>(string key, Func<Tout> function)
        {
            _manager.BeginCheckpoint(key);
            Tout output;
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
                output = Load<Tout>();
            }
            _manager.EndCheckpoint();
            return output;
        }

        public TOutput RunCheckpoint<TInput, TOutput>(
            string key,
            ISandboxedCheckpoint<TInput, TOutput> wrapper,
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
                IDirectoryLocation tempDir = GetCleanTempDirectory(key, status == CheckpointStatus.Run);
                _logger.Info("Running checkpoint {0}: {1}", checkpointNumber, checkpoint);
                Benchmark benchmark = new Benchmark();
                output = wrapper.Run(input, tempDir);
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

        public TOutput RunCheckpoint<TInput, TOutput>(
            string key,
            ISandboxedCheckpoint<TInput, TOutput> wrapper,
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
                IDirectoryLocation tempDir = GetCleanTempDirectory(key, status == CheckpointStatus.Run);
                _logger.Info("Running checkpoint {0}: {1}", checkpointNumber, checkpoint);
                Benchmark benchmark = new Benchmark();
                output = wrapper.Run(input, tempDir);
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

        public Tout RunCheckpoint<Tout>(string key, Func<IDirectoryLocation, Tout> function)
        {
            _manager.BeginCheckpoint(key);
            Tout output;
            bool hasExistingSerializedOutput = HasExistingSerializedOutput();
            CheckpointStatus status = _manager.GetCheckpointStatus(hasExistingSerializedOutput);
            if (status.ShouldRun())
            {
                string checkpoint = GetCheckpointName();
                string checkpointNumber = GetCheckpointNumber();
                IDirectoryLocation tempDir = GetCleanTempDirectory(checkpoint, status == CheckpointStatus.Run);
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
                output = Load<Tout>();
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

        private IDirectoryLocation GetCleanTempDirectory(string key, bool makeClean)
        {
            if (makeClean)
                return _placeForTempFolders.CreateEmptySubdirectory(key.RemoveWhiteSpace() + "Temp");
            else
                return _placeForTempFolders.GetDirectoryLocation(key.RemoveWhiteSpace() + "Temp");
        }

        public void Dispose()
        {
            if (_disposed) return;
            _manager.Dispose();
            _disposed = true;
        }
    }
}
