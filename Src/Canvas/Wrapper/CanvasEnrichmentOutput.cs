using System;
using Illumina.Common.FileSystem;

namespace Canvas.Wrapper
{
    public class CanvasEnrichmentOutput : ICanvasOutput
    {
        public CanvasOutput CanvasOutput { get; }
        public IFileLocation BinSize { get; }
        public IFileLocation NormalBinned { get; }
        public IFileLocation UnsmoothedCnd { get; }

        public CanvasEnrichmentOutput(CanvasOutput canvasOutput,
            IFileLocation binSize = null,
            IFileLocation normalBinned = null,
            IFileLocation unsmoothedCnd = null)
        {
            CanvasOutput = canvasOutput;
            BinSize = binSize;
            NormalBinned = normalBinned;
            UnsmoothedCnd = unsmoothedCnd;
        }

        public static CanvasEnrichmentOutput GetFromStub(IFileLocation stub, bool includeIntermediateResults)
        {
            var canvasOutput = CanvasOutput.GetFromStub(stub, includeIntermediateResults);
            if (!includeIntermediateResults)
                return new CanvasEnrichmentOutput(canvasOutput);
            IFileLocation binSize = stub.AppendName(".binsize");
            IFileLocation normalBinned = stub.AppendName(".binned");
            IFileLocation unsmoothedCnd = stub.AppendName(".unsmoothed.cnd");
            return new CanvasEnrichmentOutput(canvasOutput, binSize, normalBinned, unsmoothedCnd);
        }

        public void Move(IFileLocation fileNameStub, bool includeIntermediateResults, Action<IFileLocation, IFileLocation> move)
        {
            CanvasOutput.Move(fileNameStub, includeIntermediateResults, move);
            CanvasEnrichmentOutput destination = GetFromStub(fileNameStub, includeIntermediateResults);
            BinSize.MoveIfNotNull(destination.BinSize, move);
            NormalBinned.MoveIfNotNull(destination.NormalBinned, move);
            UnsmoothedCnd.MoveIfNotNull(destination.UnsmoothedCnd, move);
        }
    }
}