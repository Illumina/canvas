using CanvasCommon.Visualization;
using Illumina.Common;
using Isas.ClassicBioinfoTools.KentUtils;
using Isas.Framework.Logging;
using Isas.Framework.WorkManagement;
using Isas.Framework.WorkManagement.CommandBuilding;
using Isas.SequencingFiles;

namespace CanvasPedigreeCaller.Visualization
{
    public class CoverageVisualizationWriterFactory
    {
        private readonly ILogger _logger;
        private readonly IWorkDoer _workDoer;
        private readonly ICommandManager _commandManager;
        private readonly GenomeMetadata _genome;

        public CoverageVisualizationWriterFactory(ILogger logger, IWorkDoer workDoer, ICommandManager commandManager, GenomeMetadata genome)
        {
            _logger = logger;
            _workDoer = workDoer;
            _commandManager = commandManager;
            _genome = genome;
        }

        public ICoverageBigWigWriter CreateBinCoverageBigWigWriter(IBedGraphWriter bedGraphWriter)
        {
            var calculator = new NormalizedBinsCoverageCalculator();
            var bedGraphWriterFacade = new CoverageBedGraphWriter(bedGraphWriter, calculator);
            return new CoverageBigWigWriter(_logger, bedGraphWriterFacade, GetConverter(), _genome);
        }

        public CoverageBedGraphWriter CreateSegmentBedGraphWriter(IBedGraphWriter bedGraphWriter)
        {
            var calculator = new NormalizedSegmentsCoverageCalculator();
            var bedGraphWriterFacade = new CoverageBedGraphWriter(bedGraphWriter, calculator);
            return bedGraphWriterFacade;
        }


        private IBedGraphToBigWigConverter GetConverter()
        {
            if (CrossPlatform.IsThisLinux())
            {
                return new FormatConverterFactory(_logger, _workDoer, _commandManager).GetBedGraphToBigWigConverter();
            }
            return new NullBedGraphToBigWigConverter(_logger, "BedGraph to BigWig conversion unavailable on Windows.");
        }
    }
}