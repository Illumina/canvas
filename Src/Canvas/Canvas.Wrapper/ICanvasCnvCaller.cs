using Illumina.Common.FileSystem;
using Isas.Framework.DataTypes;

namespace Canvas.Wrapper
{
    public interface ICanvasCnvCaller<TCanvasInput, TCanvasOutput>
    {
        SampleSet<TCanvasOutput> Run(SampleSet<TCanvasInput> inputs, IDirectoryLocation sandbox);
    }
}