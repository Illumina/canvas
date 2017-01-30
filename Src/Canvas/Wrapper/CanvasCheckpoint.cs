using System;

namespace Illumina.SecondaryAnalysis.VariantCalling.StructuralVariants.Canvas
{
    public interface ICanvasCheckpointInput
    {
        GenomeMetadata GenomeMetadata { get; }
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
        private readonly bool _includeIntermediateResults;
        private readonly Func<IFileLocation, bool, TCanvasOutput> _getFromStub;

        public CanvasOutputNamingConventionFactory(ICanvasAnnotationFileProvider annotationFileProvider, bool includeIntermediateResults, Func<IFileLocation, bool, TCanvasOutput> getFromStub)
        {
            _annotationFileProvider = annotationFileProvider;
            _includeIntermediateResults = includeIntermediateResults;
            _getFromStub = getFromStub;
        }

        public CanvasOutputNamingConvention<TCanvasInput, TCanvasOutput> GetCanvasOutputNamingConvention(IFileLocation fileNameStub)
        {
            return new CanvasOutputNamingConvention<TCanvasInput, TCanvasOutput>(fileNameStub, _annotationFileProvider, _includeIntermediateResults, _getFromStub);
        }
    }

    public class CanvasOutputNamingConvention<TCanvasInput, TCanvasOutput> : ILoadingConvention<TCanvasInput, TCanvasOutput> where TCanvasInput : ICanvasCheckpointInput where TCanvasOutput : ICanvasOutput
    {
        private readonly IFileLocation _fileNameStub;
        private readonly ICanvasAnnotationFileProvider _annotationFileProvider;
        private readonly bool _includeIntermediateResults;
        private readonly Func<IFileLocation, bool, TCanvasOutput> _getFromStub;

        public CanvasOutputNamingConvention(IFileLocation fileNameStub, ICanvasAnnotationFileProvider annotationFileProvider, bool includeIntermediateResults, Func<IFileLocation, bool, TCanvasOutput> getFromStub)
        {
            _fileNameStub = fileNameStub;
            _annotationFileProvider = annotationFileProvider;
            _includeIntermediateResults = includeIntermediateResults;
            _getFromStub = getFromStub;
        }

        public TCanvasOutput Load(TCanvasInput resequencingInput)
        {
            if (!_annotationFileProvider.IsSupported(resequencingInput.GenomeMetadata)) return default(TCanvasOutput);
            return _getFromStub(_fileNameStub, _includeIntermediateResults);
        }

        public void Move(TCanvasOutput source, Action<IFileLocation, IFileLocation> move)
        {
            if (source == null) return;
            source.Move(_fileNameStub, _includeIntermediateResults, move);
        }
    }
}
