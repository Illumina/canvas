using Isas.ClassicBioinfoTools.KentUtils;
using Isas.Framework.Logging;
using Isas.SequencingFiles;

namespace CanvasPedigreeCaller.Visualization
{
    public class CoverageBigWigWriterFactory
    {
        private readonly ILogger _logger;
        private readonly IBedGraphToBigWigConverter _converter;
        private readonly GenomeMetadata _genome;

        public CoverageBigWigWriterFactory(ILogger logger, IBedGraphToBigWigConverter converter, GenomeMetadata genome)
        {
            _logger = logger;
            _converter = converter;
            _genome = genome;
        }

        public ICoverageBigWigWriter Create()
        {
            var calculator = new NormalizedCoverageCalculator();
            var roundingCalculator = new RoundingBedGraphCalculator(calculator, 4);
            var bedGraphWriterFacade = new BedGraphWriterFacade(roundingCalculator);
            return new CoverageBigWigWriter(_logger, bedGraphWriterFacade, _converter, _genome);
        }
    }
}