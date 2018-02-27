using System;
using Illumina.Common.FileSystem;

namespace Canvas.Wrapper
{
    public interface ICanvasOutput
    {
        void Move(IFileLocation fileNameStub, Action<IFileLocation, IFileLocation> move);
    }
}