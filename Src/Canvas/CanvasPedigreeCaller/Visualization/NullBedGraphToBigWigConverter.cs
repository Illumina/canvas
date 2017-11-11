using Illumina.Common.FileSystem;
using Isas.ClassicBioinfoTools.KentUtils;
using Isas.Framework.Logging;
using Isas.SequencingFiles;

namespace CanvasPedigreeCaller.Visualization
{
    public class NullBedGraphToBigWigConverter : IBedGraphToBigWigConverter
    {
        private readonly ILogger _logger;
        private readonly string _reasonUnavailable;

        public NullBedGraphToBigWigConverter(ILogger logger, string reasonUnavailable)
        {
            _logger = logger;
            _reasonUnavailable = reasonUnavailable;
        }

        public IFileLocation Convert(IFileLocation sourceFile, GenomeMetadata genomeMetadata, IDirectoryLocation outputDirectory)
        {
            _logger.Warn($"Not coverting {sourceFile} to BigWig. {_reasonUnavailable}");
            return null;
        }
    }
}