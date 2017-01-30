using System;

namespace Illumina.SecondaryAnalysis.VariantCalling.StructuralVariants.Canvas
{
    public interface ICanvasOutput
    {
        void Move(IFileLocation fileNameStub, bool includeIntermediateResults, Action<IFileLocation, IFileLocation> move);
    }
}