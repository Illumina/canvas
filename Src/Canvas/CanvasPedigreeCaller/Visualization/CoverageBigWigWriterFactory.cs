using Isas.ClassicBioinfoTools.KentUtils;
using Isas.Framework.Logging;
using Isas.SequencingFiles;

namespace CanvasPedigreeCaller.Visualization
{
    public class CoverageBigWigWriterFactory
    {
        private readonly ILogger _logger;
        private readonly BedGraphToBigWigConverter _converter;
        private readonly GenomeMetadata _genome;

        public CoverageBigWigWriterFactory(ILogger logger, BedGraphToBigWigConverter converter, GenomeMetadata genome)
        {
            _logger = logger;
            _converter = converter;
            _genome = genome;
        }

        public ICoverageBigWigWriter Create()
        {
            var calculator = new NormalizedCoverageCalculator();
            var bedGraphWriter = new NormalizedCoverageBedGraphWriter(calculator, 4);
            var bedGraphWriterFacade = new NormalizedCoverageWriterFacade(bedGraphWriter);
            return new CoverageBigWigWriter(_logger, bedGraphWriterFacade, _converter, _genome);
        }
    }
}