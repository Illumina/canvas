using System;
using Illumina.Common.FileSystem;
using Isas.Framework.Checkpointing;
using Isas.Framework.Checkpointing.Legacy;
using Isas.Framework.DataTypes;
using Isas.SequencingFiles;

namespace Canvas.Wrapper
{
    public interface ICanvasCheckpointInput
    {
        GenomeMetadata GenomeMetadata { get; }
    }

    public interface ILoadableResultCheckpoint<in T, in TInput, TResult> : IAbstractCheckpoint<TInput, TResult>, ISandboxedCheckpoint<TInput, TResult>
    {
        ILoadingConvention<TInput, TResult> GetLoadingConvention(T fileNameConventionSetting);
    }

    public interface ICanvasCheckpoint<TCanvasInput, TCanvasOutput> : ILoadableResultCheckpoint<SampleStubNamingConvention, SampleSet<TCanvasInput>, SampleSet<TCanvasOutput>> where TCanvasInput : ICanvasCheckpointInput where TCanvasOutput : ICanvasOutput
    {
    }

    public class CanvasCheckpoint<TCanvasInput, TCanvasOutput> : ICanvasCheckpoint<TCanvasInput, TCanvasOutput> where TCanvasInput : ICanvasCheckpointInput where TCanvasOutput : ICanvasOutput
    {
        private readonly ICanvasCnvCaller<TCanvasInput, TCanvasOutput> _canvasCnvCaller;
        private readonly CanvasOutputNamingConventionFactory<TCanvasInput, TCanvasOutput> _canvasOutputNamingConventionFactory;

        public CanvasCheckpoint(ICanvasCnvCaller<TCanvasInput, TCanvasOutput> canvasCnvCaller, CanvasOutputNamingConventionFactory<TCanvasInput, TCanvasOutput> canvasOutputNamingConventionFactory)
        {
            _canvasCnvCaller = canvasCnvCaller;
            _canvasOutputNamingConventionFactory = canvasOutputNamingConventionFactory;
        }

        public SampleSet<TCanvasOutput> Run(SampleSet<TCanvasInput> inputs, IDirectoryLocation sandbox)
        {
            return _canvasCnvCaller.Run(inputs, sandbox);
        }

        public string GetStepName()
        {
            return "Detect CNV";
        }

        public ILoadingConvention<SampleSet<TCanvasInput>, SampleSet<TCanvasOutput>> GetLoadingConvention(SampleStubNamingConvention fileNameConventionSetting)
        {
            return new SampleSetLoadingConvention<TCanvasInput, TCanvasOutput>(fileNameConventionSetting, (fileNameStub) => _canvasOutputNamingConventionFactory.GetCanvasOutputNamingConvention(fileNameStub));
        }
    }

    public class CanvasOutputNamingConventionFactory<TCanvasInput, TCanvasOutput> where TCanvasInput : ICanvasCheckpointInput where TCanvasOutput : ICanvasOutput
    {
        private readonly ICanvasAnnotationFileProvider _annotationFileProvider;
        private readonly Func<IFileLocation, TCanvasOutput> _getFromStub;

        public CanvasOutputNamingConventionFactory(ICanvasAnnotationFileProvider annotationFileProvider, Func<IFileLocation, TCanvasOutput> getFromStub)
        {
            _annotationFileProvider = annotationFileProvider;
            _getFromStub = getFromStub;
        }

        public CanvasOutputNamingConvention<TCanvasInput, TCanvasOutput> GetCanvasOutputNamingConvention(IFileLocation fileNameStub)
        {
            return new CanvasOutputNamingConvention<TCanvasInput, TCanvasOutput>(fileNameStub, _annotationFileProvider, _getFromStub);
        }
    }

    public class CanvasOutputNamingConvention<TCanvasInput, TCanvasOutput> : ILoadingConvention<TCanvasInput, TCanvasOutput> where TCanvasInput : ICanvasCheckpointInput where TCanvasOutput : ICanvasOutput
    {
        private readonly IFileLocation _fileNameStub;
        private readonly ICanvasAnnotationFileProvider _annotationFileProvider;
        private readonly Func<IFileLocation, TCanvasOutput> _getFromStub;

        public CanvasOutputNamingConvention(IFileLocation fileNameStub, ICanvasAnnotationFileProvider annotationFileProvider, Func<IFileLocation, TCanvasOutput> getFromStub)
        {
            _fileNameStub = fileNameStub;
            _annotationFileProvider = annotationFileProvider;
            _getFromStub = getFromStub;
        }

        public TCanvasOutput Load(TCanvasInput resequencingInput)
        {
            if (!_annotationFileProvider.IsSupported(resequencingInput.GenomeMetadata)) return default(TCanvasOutput);
            return _getFromStub(_fileNameStub);
        }

        public void Move(TCanvasOutput source, Action<IFileLocation, IFileLocation> move)
        {
            if (source == null) return;
            source.Move(_fileNameStub, move);
        }
    }

    public class SampleSetLoadingConvention<TIn, T> : ILoadingConvention<SampleSet<TIn>, SampleSet<T>>
    {
        private readonly NamingConventionGenerator _loadingConventionGenerator;
        private readonly SampleStubNamingConvention _convention;

        public delegate ILoadingConvention<TIn, T> NamingConventionGenerator(IFileLocation fileNameStub);

        public SampleSetLoadingConvention(SampleStubNamingConvention sampleStubNamingConvention, NamingConventionGenerator loadingConventionGenerator)
        {
            _loadingConventionGenerator = loadingConventionGenerator;
            _convention = sampleStubNamingConvention;
        }

        public SampleSet<T> Load(SampleSet<TIn> input)
        {
            var set = new SampleSet<T>();
            foreach (var sample in input)
            {
                IFileLocation stub = _convention(sample.Key);
                ILoadingConvention<TIn, T> namingConvention = _loadingConventionGenerator(stub);
                T output = namingConvention.Load(sample.Value);
                set.Add(sample.Key, output);
            }
            return set;
        }

        public void Move(SampleSet<T> source, Action<IFileLocation, IFileLocation> move)
        {
            foreach (var sampleInfo in source.Keys)
            {
                IFileLocation stub = _convention(sampleInfo);
                ILoadingConvention<TIn, T> namingConvention = _loadingConventionGenerator(stub);
                namingConvention.Move(source[sampleInfo], move);
            }
        }
    }

    public static class LoadingAndMovingCheckpointer
    {
        public static TResult RunCheckpoint<TInput, TResult>(this ICheckpointRunner runner, string name, Func<TInput, IDirectoryLocation, TResult> run, TInput input, ILoadingConvention<TInput, TResult> loadingConvention)
        {
            return runner.RunCheckpoint(name, (dir, mover) =>
            {
                var result = run(input, dir);
                loadingConvention.Move(result, mover.Move);
                return loadingConvention.Load(input);
            }, () => loadingConvention.Load(input));
        }

        public static TResult RunCheckpoint<TResult>(this ICheckpointRunner runner, string name, Func<TResult> run, INamingConvention<TResult> namingConvention)
        {
            return runner.RunCheckpoint(name, (dir, mover) =>
            {
                var result = run();
                return namingConvention.Move(result, mover.Move);
            });
        }

        public static TResult RunCheckpoint<TResult>(this ICheckpointRunner runner, string name, Func<IDirectoryLocation, TResult> run, INamingConvention<TResult> namingConvention)
        {
            return runner.RunCheckpoint(name, (dir, mover) =>
            {
                var result = run(dir);
                return namingConvention.Move(result, mover.Move);
            });
        }

        public static TOutput RunCheckpoint<TConvention, TInput, TOutput>(this ICheckpointRunner runner,
    ILoadableResultCheckpoint<TConvention, TInput, TOutput> wrapper,
    TInput input,
    TConvention convention)
        {
            return runner.RunCheckpoint(wrapper, input, wrapper.GetLoadingConvention(convention));
        }

        public static TOutput RunCheckpoint<TConvention, TInput, TOutput>(this ICheckpointRunner runner,
    IMoveableResultCheckpoint<TConvention, TInput, TOutput> wrapper,
    TInput input,
    TConvention convention)
        {
            return runner.RunCheckpoint(wrapper, input, wrapper.GetNamingConvention(convention));
        }

        public static TOutput RunCheckpoint<TInput, TOutput>(this ICheckpointRunner runner,
    IAbstractCheckpoint<TInput, TOutput> wrapper,
    TInput input,
    ILoadingConvention<TInput, TOutput> loadingConvention)
        {
            return runner.RunCheckpoint(wrapper.GetStepName(), wrapper.Run, input, loadingConvention);
        }

        public static TOutput RunCheckpoint<TInput, TOutput>(this ICheckpointRunner runner,
IAbstractCheckpoint<TInput, TOutput> wrapper,
TInput input,
INamingConvention<TOutput> namingConvention)
        {
            return runner.RunCheckpoint(wrapper.GetStepName(), dir => wrapper.Run(input, dir), namingConvention);
        }
    }
}
