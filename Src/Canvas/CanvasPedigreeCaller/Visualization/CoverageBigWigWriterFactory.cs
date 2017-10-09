using Isas.ClassicBioinfoTools.KentUtils;
using Isas.SequencingFiles;

namespace CanvasPedigreeCaller.Visualization
{
    public class CoverageBigWigWriterFactory
    {
        private readonly BedGraphToBigWigConverter _converter;
        private readonly GenomeMetadata _genome;

        public CoverageBigWigWriterFactory(BedGraphToBigWigConverter converter, GenomeMetadata genome)
        {
            _converter = converter;
            _genome = genome;
        }

        public ICoverageBigWigWriter Create()
        {
            var calculator = new NormalizedCoverageCalculator();
            var bedGraphWriter = new NormalizedCoverageBedGraphWriter(calculator, 4);
            var bedGraphWriterFacade = new NormalizedCoverageWriterFacade(bedGraphWriter);
            return new CoverageBigWigWriter(bedGraphWriterFacade, _converter, _genome);
        }
    }
}