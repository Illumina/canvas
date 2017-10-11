using System.Collections.Generic;
using CanvasCommon;
using Illumina.Common.FileSystem;
using Isas.ClassicBioinfoTools.KentUtils;
using Isas.Framework.Logging;
using Isas.Framework.Utilities;
using Isas.SequencingFiles;

namespace CanvasPedigreeCaller.Visualization
{
    public class CoverageBigWigWriter : ICoverageBigWigWriter
    {
        private readonly ILogger _logger;
        private readonly ICoverageBedGraphWriter _writer;
        private readonly BedGraphToBigWigConverter _converter;
        private readonly GenomeMetadata _genome;

        public CoverageBigWigWriter(ILogger logger, ICoverageBedGraphWriter writer, BedGraphToBigWigConverter converter, GenomeMetadata genome)
        {
            _logger = logger;
            _writer = writer;
            _converter = converter;
            _genome = genome;
        }

        public IFileLocation Write(IReadOnlyList<CanvasSegment> segments, IDirectoryLocation output)
        {
            _logger.Info($"Begin writing bedgraph file at '{output}'");
            var benchmark = new Benchmark();
            var bedGraph = output.GetFileLocation("coverage.bedgraph");
            _writer.Write(segments, bedGraph);
            _logger.Info($"Finished writing bedgraph file at '{output}'. Elapsed time: {benchmark.GetElapsedTime()}");
            _logger.Info($"Begin conversion of '{output}' to bigwig file");
            benchmark = new Benchmark();
            var bigWigConverterOutput = output.CreateSubdirectory("BigWigConverter");
            var bigwigFile = _converter.Convert(bedGraph, _genome, bigWigConverterOutput);
            _logger.Info($"Finished conversion from bedgraph file at '{output}' to bigwig file at '{bigwigFile}'. Elapsed time: {benchmark.GetElapsedTime()}");
            return bigwigFile;
        }
    }
}