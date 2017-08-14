using System.Threading.Tasks;
using Isas.Framework.Checkpointing;
using Isas.Framework.Checkpointing.Legacy;
using Isas.Framework.DataTypes;

namespace Canvas.Wrapper
{
    public interface ICanvasWorker<TCanvasInput, TCanvasOutput> where TCanvasOutput : ICanvasOutput
    {
        SampleSet<TCanvasOutput> Run(SampleSet<TCanvasInput> inputs, SampleStubNamingConvention sampleStubNamingConvention, ICheckpointRunner checkpointRunner);
        Task<SampleSet<TCanvasOutput>> RunAsync(SampleSet<TCanvasInput> inputs, SampleStubNamingConvention sampleStubNamingConvention, ICheckpointRunner checkpointRunner);
    }

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
}
